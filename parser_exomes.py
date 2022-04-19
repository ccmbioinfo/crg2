import argparse
import ast
import glob
import os
import sys
import pandas as pd
import subprocess
from collections import namedtuple


"""
Usage: python parser.py -f <input_sample.tsv> -d <absolute path to create directories> -s [mapped|recal|fastq|decoy_rm] -b <bioinformaticians>
Parses analyses request csv from Stager and sets up necessary directories (under 2nd argument), files as below:
1. create family and directory passed as "-s"
2. symlink BAM files if "-s" is not fastq
3. copy config.yaml, pbs_config.yaml and dnaseq_cluster.pbs from crg2 repo and replace necessary string
4. searches hpf for input files uploaded via MinIO, or from previously run analyses
5. create units.tsv and samples.tsv for snakemake
6. submit job if all the above goes well
"""

crg2_dir = os.path.dirname(os.path.realpath(sys.argv[0]))


def input_file(input_path):
    """Given hpf path, find input files"""
    fq1 = sorted(glob.glob(os.path.join(input_path, "*_R1*.fastq.gz"))) + sorted(
        glob.glob(os.path.join(input_path, "*_1.fastq.gz"))
    )
    fq2 = sorted(glob.glob(os.path.join(input_path, "*_R2*.fastq.gz"))) + sorted(
        glob.glob(os.path.join(input_path, "*_2.fastq.gz"))
    )
    bam = glob.glob(os.path.join(input_path, "*bam"))
    cram = glob.glob(os.path.join(input_path, "*cram"))
    # prioritize fastq as input, then bam, then cram
    input_types = ("fastq", "bam", "cram")
    input = dict.fromkeys(input_types, "")
    if len(fq1) != 0:
        if len(fq1) == len(fq2):
            input["fastq"] = {"R1": fq1, "R2": fq2}
        else:
            print(f"Number of R1 and R2 files in {input_path} does not match")
            input = None
    elif len(bam) != 0:
        if len(bam) == 1:
            input["bam"] = bam[0]
        else:
            print(f"Multiple bams found under {input_path}")
            input = None
    elif len(cram) != 0:
        if len(cram) == 1:
            input["cram"] = cram[0]
        else:
            print(f"Multiple crams found under {input_path}")
            input = None
    else:
        input = None

    return input


def check_minio_dccforge(family, participant, sequencing_id):
    """Given family and participant codenames, find input files on the hpf"""
    DCCFORGE_UPLOADS = "/hpf/largeprojects/ccm_dccforge/dccforge/uploads/*/"
    MINIO_MIRROR = (
        "/hpf/largeprojects/ccm_dccforge/dccforge/uploads/MinIO_G4RD_Mirror/*/"
    )
    DPLM = "/hpf/largeprojects/ccmbio_ephemeral/C4R/"
    dccforge_input = input_file(f"{DCCFORGE_UPLOADS}/{family}_{participant}*/")
    minio_input = input_file(f"{MINIO_MIRROR}/{family}_{participant}*/")
    # some datasets are transferred to the ccmbio_ephemeral space by DPLM
    dplm_input = None
    if sequencing_id != None:
        dplm_input = input_file(f"{DPLM}/*{sequencing_id}*/")
    if not dccforge_input and not minio_input and not dplm_input:
        input = None
    else:
        if minio_input:
            input = minio_input
        elif dplm_input:
            input = dplm_input
        else:
            input = dccforge_input
    return input


def find_input(family, participant, sequencing_id):
    """Given family and participant codenames, check if individuals have already been analyzed"""
    RESULTS_DIR = "/hpf/largeprojects/ccm_dccforge/dccforge/results/*/"
    # check uploads folders first, as participant may have been resequenced
    input = check_minio_dccforge(family, participant, sequencing_id)
    if not input:
        # if not in uploads, this is a re-analysis using bams from results directory
        bam = glob.glob(f"{RESULTS_DIR}/{family}/{family}_{participant}.bam")
        if len(bam) == 1:
            input_types = ("fastq", "bam", "cram")
            input = dict.fromkeys(input_types, "")
            input["bam"] = bam
        # if not in results, may be missing
        else:
            input = None

    return input


def valid_dir(dir):
    """Check that directory provided exists"""
    if os.path.isdir(dir):
        return dir
    else:
        message = f"{dir} path does not exist. Please provide absolute path"
        print(message)
        raise argparse.ArgumentTypeError(message)


def valid_file(filename):
    """Check that file exists and is not empty"""
    if not os.path.isfile(filename):
        message = f"{filename} file does not exist"
        raise argparse.ArgumentTypeError(message)
    else:
        if not os.path.getsize(filename) > 0:
            message = f"{filename} file is empty"
            print(message)
            raise argparse.ArgumentTypeError(message)
    return filename


def assign_exomes(requested, bioinfos):
    """Assign analyses to bioinfos"""
    num_analyses = len(requested)
    remainder = num_analyses % len(bioinfos)
    num_analyses_per = (num_analyses - remainder) / len(bioinfos)

    assignees = []
    for b in bioinfos:
        assignees.append([b] * (int(num_analyses_per)))
    assignees.append([bioinfos[0]] * remainder)
    # flatten assignees sublists in assignees list
    assignees = [item for sublist in assignees for item in sublist]
    requested["assignee"] = assignees

    return requested


def setup_directories(family, sample_list, filepath, step=None):
    d = os.path.join(filepath, family)
    if os.path.isdir(d):
        print(f"{family} analysis directory already exists")
        inputs_flag = False
    else:
        os.mkdir(d)
        # copy config.yaml, pbs_config.yaml, and dnaseq_cluster.pbs
        for i in ["config.yaml", "pbs_profile/pbs_config.yaml", "dnaseq_cluster.pbs"]:
            cmd = ["cp", os.path.join(crg2_dir, i), d]
            subprocess.check_call(cmd)

        # replace family ID in config.yaml & dnaseq_cluster.pbs
        replace = f"s/NA12878/{family}/"
        config = os.path.join(d, "config.yaml")
        if os.path.isfile(config):
            cmd = ["sed", "-i", replace, config]
            subprocess.check_call(cmd)
        replace = f"s/crg2_pbs/{family}/"
        pbs = os.path.join(d, "dnaseq_cluster.pbs")
        if os.path.isfile(pbs):
            cmd = ["sed", "-i", replace, pbs]
            subprocess.check_call(cmd)

        # check to see if each sample is associated with input file
        inputs_flag = check_inputs(sample_list)

    return inputs_flag


def check_inputs(sample_list):
    """Checks if each sample (a named tuple) has an input file (fastq, bam, or cram)"""
    check_inputs = True
    for i in sample_list:
        if not i.fq1 and not i.fq2 and not i.bam and not i.cram:
            check_inputs = False
    return check_inputs


def input_types(sample):
    """Check if input type is fastq or bam"""
    if sample.fq1:
        return "fastq"
    elif sample.cram:
        return "cram"
    else:
        return "bam"


def write_proj_files(sample_list, filepath):
    units, samples = os.path.join(filepath, "units.tsv"), os.path.join(
        filepath, "samples.tsv"
    )
    with open(units, "w") as u, open(samples, "w") as s:
        s.writelines("sample\n")
        u.writelines("sample\tplatform\tfq1\tfq2\tbam\tcram\n")
        for i in sample_list:
            s.writelines(f"{i.sample}\n")
            u.writelines(
                f"{i.sample}\t{i.platform}\t{i.fq1}\t{i.fq2}\t{i.bam}\t{i.cram}\n"
            )


def submit_jobs(directory, family):
    if not os.path.isdir(directory):
        print(f"{directory} not found. Exiting!")
        exit()
    pbs = "dnaseq_cluster.pbs"
    if os.path.isfile(os.path.join(directory, pbs)):
        cwd = os.getcwd()
        os.chdir(directory)
        cmd = ["qsub", pbs]
        jobid = subprocess.check_output(cmd).decode("UTF-8").rstrip()
        print(f"Submitted snakemake job for {family}: {jobid}")
        os.chdir(cwd)
    else:
        print(f"{pbs} not found, not submitting jobs for {family}.")


if __name__ == "__main__":
    description = "Reads sample info from csv file (-f) and creates directory (-d) necessary to start crg2 from step (-s) requested."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-f",
        "--file",
        type=valid_file,
        required=True,
        help="Analyses csv output from STAGER",
    )

    parser.add_argument(
        "-d",
        "--dir",
        type=valid_dir,
        required=True,
        metavar="path",
        help="Absolute path where crg2 directory structure will be created under familyID as base directory",
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
        "-b",
        "--bioinfos",
        nargs="+",
        type=str,
        required=True,
        help="Names of bioinformaticians who will be assigned cases separated by spaces, e.g. MC NH PX",
    )

    args = parser.parse_args()
    requested = pd.read_csv(args.file)
    inputs_list = []

    Units = namedtuple("Units", "sample platform fq1 fq2 bam cram")
    platform = "ILLUMINA"
    projects = {}
    for index, row in requested.iterrows():
        family = list(set(ast.literal_eval(row["family_codenames"])))
        # Will need to modify for multi-family analyses (i.e. cohort)
        family = family[0]
        if not family in projects:
            projects[family] = []
        family_inputs = []
        for participant, sequencing_id in zip(
            ast.literal_eval(row["participant_codenames"]),
            ast.literal_eval(row["sequencing_id"]),
        ):
            # look for input files in MinIO mirror on hpf, or in results directory
            input = find_input(family, participant, sequencing_id)
            if input:
                if input["fastq"]:
                    fq1, fq2 = (",").join(input["fastq"]["R1"]), (",").join(
                        input["fastq"]["R2"]
                    )
                else:
                    fq1, fq2 = "", ""
                bam, cram = input["bam"], input["cram"]
                projects[family].append(
                    Units(participant, platform, fq1, fq2, bam, cram)
                )
            else:
                fq1, fq2, bam, cram = "", "", "", ""
                projects[family].append(
                    Units(participant, platform, fq1, fq2, bam, cram)
                )
                print(f"Input files are missing for {family}_{participant}")

            family_inputs.append(input)
        inputs_list.append(family_inputs)

    for i in projects:

        sample_list = projects[i]
        filepath = os.path.join(args.dir, i)

        # create project directory
        # copy config.yaml, dnaseq_cluster.pbs & replace family id; copy pbs_config.yaml
        print(f"\nStarting to setup project directories for family: {i}")
        submit_flag = setup_directories(i, sample_list, args.dir, args.step)

        write_proj_files(sample_list, filepath)
        if submit_flag:
            print(f"\nSubmitting job for family: {i}")
            # submit_jobs(filepath, i)
        else:
            print(f"Job for {i} not submitted")

    requested["input"] = inputs_list
    requested = assign_exomes(requested, args.bioinfos)
    requested_name = args.file.strip(".csv")
    requested.to_csv(f"{requested_name}_with_input_assigned.csv", index=False)
