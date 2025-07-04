import sys, os, subprocess, glob
from collections import namedtuple
import argparse
import pandas as pd

"""
Usage: python parser.py -f <fastq_input_sample.tsv> -a <stager_analyses_report.csv> -d <absolute path to create directories> -s [mapped|recal|fastq|decoy_rm] -p <project ID for sequence batch, e.g. ACK23274>
Parse five-column(family,sample,fq1,fq2,bam) TSV file (1st argument) and sets up necessary directories (under 2nd argument), files as below:
1. create family and directory passed as "-s"
2. symlink BAM files if "-s" is not fastq
3. copy config_cheo_ri.yaml, slurm-config.yaml and dnaseq_slurm_cheo_ri.sh from crg2 repo and replace necessary string
4. create units.tsv and samples.tsv for snakemake
5. generate CNV report from TCAG CNV tsvs
6. submit job if all the above goes well
"""

crg2_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
folder_file_suffix = {
    "mapped": "sorted.bam",
    "recal": "bam",
    "decoy_rm": "no_decoy_reads.bam",
}


def create_symlink(src, dest):
    if not os.path.isfile(dest):
        cmd = ["ln", "-s", src, dest]
        subprocess.check_call(cmd)


def setup_directories(family, sample_list, filepath, step):
    d = filepath
    if not os.path.isdir(d):
        print(d)
        os.makedirs(d, exist_ok=True)

    # copy config_cheo_ri.yaml, slurm_config.yaml, and dnaseq_slurm_cheo_ri.sh
    for i in [
        "config_cheo_ri.yaml",
        "slurm_profile/slurm-config.yaml",
        "dnaseq_slurm_cheo_ri.sh",
    ]:
        cmd = ["cp", os.path.join(crg2_dir, i), d]
        subprocess.check_call(cmd)

    # replace family ID and pipeline in config_cheo_ri.yaml & dnaseq_slurm_cheo_ri.sh
    replace = "s/NA12878/{}/".format(family)
    pipeline = "s/wes/wgs/"
    PT_credentials = 's+PT_credentials: ""+PT_credentials: {}+'.format(
        "/srv/shared/crg2/credentials.csv"
    )
    config = os.path.join(d, "config_cheo_ri.yaml")
    if os.path.isfile(config):
        cmd = ["sed", "-i", replace, config]
        subprocess.check_call(cmd)
        cmd = ["sed", "-i", pipeline, config]
        subprocess.check_call(cmd)
        cmd = ["sed", "-i", PT_credentials, config]
        subprocess.check_call(cmd)
    replace = "s/job-name=crg2/job-name={}/".format(family)
    config_path = "s+/srv/shared/crg2/config_cheo_ri.yaml+{}/config_cheo_ri.yaml+".format(d)
    slurm = os.path.join(d, "dnaseq_slurm_cheo_ri.sh")
    if os.path.isfile(slurm):
        cmd = ["sed", "-i", replace, slurm]
        subprocess.check_call(cmd)
        cmd = ["sed", "-i", config_path, slurm]
        subprocess.check_call(cmd)

    # glob hpo
    hpo_path = os.path.expanduser("/srv/shared/metadata/HPO")
    try:
        hpo = glob.glob(("{}/{}_HPO*.txt").format(hpo_path, family))[-1] # get most recent HPO file
    except IndexError:
        print(f"No HPO files found for family: {family}")
        hpo = ""

    if os.path.isfile(config):
        replace = 's+hpo: ""+hpo: "{}"+'.format(hpo)
        cmd = ["sed", "-i", replace, config]
        subprocess.check_call(cmd)

    # glob ped
    ped_path = os.path.expanduser("/srv/shared/metadata/pedigrees")
    try:
        ped = glob.glob(("{}/{}*ped").format(ped_path, family))[-1] # get most recent ped file
    except IndexError:
        print(f"No ped files found for family: {family}")
        ped = ""

    if os.path.isfile(config):
        pedi = pd.read_csv(
            ped,
            sep=" ",
            header=None,
            names=["fam_id", "individual_id", "pat_id", "mat_id", "sex", "phenotype"],
        )
        for index, row in pedi.iterrows():
            individual_id = str(row.individual_id)
            pat_id = str(row.pat_id)
            mat_id = str(row.mat_id)
            # individual_id, pat_id and mat_id cannot be both numeric and single number
            if not (
                (individual_id.isnumeric() and len(individual_id) == 1)
                | (pat_id.isnumeric() and len(pat_id) == 1)
                | (mat_id.isnumeric() and len(mat_id) == 1)
            ):
                # write and check sample
                write_sample(filepath, family)
                samples = os.path.join(filepath, "samples.tsv")
                samples = pd.read_csv(samples, dtype=str)
                samples = list(samples.iloc[:, 0])
                samples = [family + "_" +  s for s in samples]
                samples.sort()
                print(f"samples to be processed: {samples}")
                ped_samples = [individual_id, pat_id, mat_id]
                ped_samples.sort()
                if all(s in samples for s in ped_samples):
                    replace = 's+ped: ""+ped: "{}"+'.format(ped)
                    cmd = ["sed", "-i", replace, config]
                    subprocess.check_call(cmd)
                    break
                else:
                    print(
                        f"Individuals in ped file do not match with samples.tsv, double check: {ped}."
                    )
                    return False
            else:
                print(
                    f"Ped file is either not a trio or not linked, double check: {ped}."
                )
                return False

    # bam: start after folder creation and symlink for step
    if step in ["mapped", "decoy_rm", "recal"]:
        if len(sample_list) == len([i for i in sample_list if i.bam]):
            start_folder = os.path.join(d, step)
            if not os.path.isdir(start_folder):
                os.mkdir(start_folder)
                for i in sample_list:
                    dest = "{}_{}.{}".format(family, i.sample, folder_file_suffix[step])
                    create_symlink(i.bam, os.path.join(start_folder, dest))
                    if os.path.isfile(i.bam + ".bai"):
                        src = i.bam + ".bai"
                        dest = dest + ".bai"
                        create_symlink(src, os.path.join(start_folder, dest))
            else:
                print(
                    f"{start_folder} directory already exists! Won't rewrite/submit jobs. "
                )
                return False
        return True

    # fastq: no directory creations required; inputs can be fastq or bam
    # set input: bam or fastq(default) in config_cheo_ri.yaml
    if step == "fastq":
        if len(sample_list) == len([i for i in sample_list if i.bam]):
            return True
        if len(sample_list) == len([i for i in sample_list if i.fq1]):
            return True


def write_sample(filepath, family):
    samples = os.path.join(filepath, "samples.tsv")
    with open(samples, "w") as s:
        s.writelines("sample\n")
        for i in sample_list:
            s.writelines(f"{i.sample}\n")


def write_units(sample_list, filepath):
    units = os.path.join(filepath, "units.tsv")
    with open(units, "w") as u:
        u.writelines("sample\tplatform\tfq1\tfq2\tbam\n")
        for i in sample_list:
            u.writelines(f"{i.sample}\t{i.platform}\t{i.fq1}\t{i.fq2}\t{i.bam}\n")


def submit_jobs(directory, family):
    if not os.path.isdir(directory):
        print(f"{directory} not found. Exiting!")
        exit()
    slurm = "dnaseq_slurm_cheo_ri.sh"
    if os.path.isfile(os.path.join(directory, slurm)):
        cwd = os.getcwd()
        os.chdir(directory)
        cmd = ["sbatch", slurm]
        jobid = subprocess.check_output(cmd).decode("UTF-8").rstrip()
        print(f"Submitted snakemake job for {family}: {jobid}")
        os.chdir(cwd)
    else:
        print(f"{slurm} not found, not submitting jobs for {family}.")


def valid_dir(dir):
    if os.path.isdir(dir):
        return dir
    else:
        message = "{} path does not exist. Please provide absolute path".format(dir)
        print(message)
        raise argparse.ArgumentTypeError(message)


def valid_file(filename):
    if not os.path.isfile(filename):
        message = "{} file does not exist".format(filename)
        print(message)
        raise argparse.ArgumentTypeError(message)
    else:
        if not os.path.getsize(filename) > 0:
            message = "{} file is empty".format(filename)
            print(message)
            raise argparse.ArgumentTypeError(message)
    return filename


def cp_cnv(project, family, filepath):
    cnvs = glob.glob(
        f"/srv/shared/data/TCAG_SRWGS_raw/{project}_cnv_*/FAM_{family}/cnv/*tsv"
    )
    os.makedirs(f"{filepath}/report/cnv", exist_ok=True)
    d = f"{filepath}/report/cnv"
    for c in cnvs:
        cmd = ["cp", c, d]
        subprocess.check_call(cmd)
    crg_dir = crg2_dir.replace("/crg2/", "/crg/")
    cmd = ["sbatch", f"{crg_dir}/merge.cnv.reports.sh", family]
    os.chdir(d)
    subprocess.check_call(cmd)


if __name__ == "__main__":
    description = """Reads sample info from TSV file (-f) and creates directory (-d) necessary to start
    crg2 from step (-s) requested.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-f",
        "--fastq_file",
        type=valid_file,
        required=True,
        help="Five column TAB-seperated file containing fastq paths for each sample sequenced in TCAG batch",
    )
    parser.add_argument(
        "-a",
        "--analyses_csv",
        type=valid_file,
        required=True,
        help="Analyses report CSV downloaded from Stager containing analyses corresponding to samples in fastq_file",
    )
    parser.add_argument(
        "-s",
        "--step",
        type=str,
        required=True,
        choices=["fastq", "mapped", "recal", "decoy_rm"],
        help="start running from this folder creation(step)",
    )
    parser.add_argument(
        "-d",
        "--dir",
        type=valid_dir,
        required=True,
        metavar="path",
        help="Absolute path where crg2 directory struture will be created under familyID as base directory",
    )
    parser.add_argument(
        "-p",
        "--project",
        type=str,
        required=True,
        help="Project ID for sequence batch, e.g. ACK23274",
    )

    projects = {}
    Units = namedtuple("Units", "sample platform fq1 fq2 bam")
    platform = "ILLUMINA"
    args = parser.parse_args()
    with open(args.fastq_file) as f:
        for i in f:
            if not i.startswith("familyID"):
                family, sample, fq1, fq2, bam = i.strip("\n").split("\t")
                if not family in projects:
                    projects[family] = []
                projects[family].append(Units(sample, platform, fq1, fq2, bam))

    analyses_report = pd.read_csv(args.analyses_csv)
    analyses_report["family_codenames"] = analyses_report["family_codenames"].apply(lambda x: x.strip("[]\'").replace("'", "").split(", ")[0])
    project_analysis =  dict(zip(analyses_report['family_codenames'], analyses_report['analysis_id']))

    for i in projects:
        sample_list = projects[i]
        analysis_id = project_analysis[i]
        filepath = os.path.join(args.dir, i, str(analysis_id))
        print(i, filepath)

        # create project directory
        # copy config_cheo_ri.yaml, dnaseq_slurm_cheo_ri.sh & replace family id; copy slurm_config.yaml
        print(f"\nStarting to setup project directories for family: {i}")
        submit_flag = setup_directories(i, sample_list, filepath, args.step)
        write_units(sample_list, filepath)
        cp_cnv(args.project, i, filepath)

        if submit_flag:
            print("submit jobs")
            submit_jobs(filepath, i)
        else:
            write_sample(filepath, i)

    print("DONE")
