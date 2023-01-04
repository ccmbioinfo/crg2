import sys, os, subprocess, glob
from collections import namedtuple
import argparse

<<<<<<< HEAD
'''
=======
"""

>>>>>>> 3c5d678d0fee2c2492072798a04fb2663bc1e247
Usage: python parser.py -f <input_sample.tsv> -d <absolute path to create directories> -s [mapped|recal|fastq|decoy_rm]
Parse five-column(family,sample,fq1,fq2,bam) TSV file (1st argument) and sets up necessary directories (under 2nd argument), files as below:
1. create family and directory passed as "-s"
2. symlink BAM files if "-s" is not fastq
3. copy config.yaml, pbs_config.yaml and dnaseq_cluster.pbs from crg2 repo and replace necessary string
4. create units.tsv and samples.tsv for snakemake
5. submit job if all the above goes well
'''

crg2_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
folder_file_suffix = {"mapped": "sorted.bam", "recal": "bam", "decoy_rm": "no_decoy_reads.bam"}
def create_symlink(src, dest):
    if not os.path.isfile(dest):
        cmd = ["ln", "-s", src, dest]
        subprocess.check_call(cmd)

def setup_directories(family, sample_list, filepath, step):
    d = os.path.join(filepath, family)
    if not os.path.isdir(d):
        os.mkdir(d)

    #copy config.yaml, pbs_config.yaml, and dnaseq_cluster.pbs
    for i in ["config.yaml", "pbs_profile/pbs_config.yaml", "dnaseq_cluster.pbs"]:
        cmd = ["cp", os.path.join(crg2_dir, i), d]
        subprocess.check_call(cmd)

    #replace family ID and pipeline in config.yaml & dnaseq_cluster.pbs
    replace = 's/NA12878/{}/'.format(family)
    pipeline = 's/wes/wgs/'
    config = os.path.join(d,"config.yaml")
    if os.path.isfile(config):
        cmd = ["sed", "-i", replace, config]
        subprocess.check_call(cmd)
        cmd = ["sed", "-i", pipeline, config]
        subprocess.check_call(cmd)
    replace = 's/crg2_pbs/{}/'.format(family)
    pbs = os.path.join(d,"dnaseq_cluster.pbs")
    if os.path.isfile(pbs):
        cmd = ["sed", "-i", replace, pbs]
        subprocess.check_call(cmd)

    #glob hpo
    hpo_path = os.path.expanduser('~/gene_data/HPO')
    hpo = glob.glob(('{}/{}_HPO.txt').format(hpo_path, family))
    if len(hpo) > 1:
        print(f"Multiple HPO files found: {hpo}. Exiting!")
        exit()
    if len(hpo) == 1 and os.path.isfile(config):
        hpo = hpo[0]
        replace = 's+hpo: ""+hpo: "{}"+'.format(hpo)
        cmd = ["sed", "-i", replace, config]
        subprocess.check_call(cmd)

    #glob ped
    ped_path = os.path.expanduser('~/gene_data/pedigrees')
    ped = glob.glob(('{}/{}*ped').format(ped_path, family))
    if len(ped) > 1:
        print(f"Multiple ped files found: {ped}. Exiting!")
        exit()
    if len(ped) == 1 and os.path.isfile(config):
        ped = ped[0]
        replace = 's+ped: ""+ped: "{}"+'.format(ped)
        cmd = ["sed", "-i", replace, config]
        subprocess.check_call(cmd)

    # bam: start after folder creation and symlink for step
    if step in ["mapped", "decoy_rm", "recal"]:
        if len(sample_list) == len([i for i in sample_list if i.bam]):

            start_folder = os.path.join(d, step)
            if not os.path.isdir(start_folder):
                os.mkdir(start_folder)
                for i in sample_list:
                    dest = "{}_{}.{}".format(family,i.sample,folder_file_suffix[step])
                    create_symlink(i.bam, os.path.join(start_folder, dest))
                    if os.path.isfile(i.bam + ".bai"):
                        src = i.bam + ".bai"
                        dest = dest + ".bai"
                        create_symlink(src, os.path.join(start_folder, dest))
            else:
                print(f"{start_folder} directory already exists! Won't rewrite/submit jobs. ")
                return False
        return True
    
    #fastq: no directory creations required; inputs can be fastq or bam
    #set input: bam or fastq(default) in config.yaml
    if step == "fastq":
        if len(sample_list) == len([i for i in sample_list if i.bam ]):
            replace = 's/input: "fastq"/input: "bam"/g'
            config = os.path.join(d,"config.yaml")
            if os.path.isfile(config):
                cmd = ["sed", "-i", replace, config]
                subprocess.check_call(cmd)
            return True
        if len(sample_list) == len([i for i in sample_list if i.fq1 ]): 
            return True
    

def write_proj_files(sample_list, filepath):
    units, samples = os.path.join(filepath, "units.tsv"), os.path.join(filepath, "samples.tsv")
    with open(units, "w") as u, open(samples, "w") as s:
        s.writelines("sample\n")
        u.writelines("sample\tplatform\tfq1\tfq2\tbam\n")
        for i in sample_list:
            s.writelines(f"{i.sample}\n")
            u.writelines(f"{i.sample}\t{i.platform}\t{i.fq1}\t{i.fq2}\t{i.bam}\n")
            

def submit_jobs(directory, family):
    if not os.path.isdir(directory):
        print(f"{directory} not found. Exiting!")
        exit()
    pbs =  "dnaseq_cluster.pbs"
    if os.path.isfile(os.path.join(directory,pbs)):
        cwd = os.getcwd()
        os.chdir(directory)
        cmd = ["qsub", pbs]
        jobid = subprocess.check_output(cmd).decode("UTF-8").rstrip()
        print(f"Submitted snakemake job for {family}: {jobid}")
        os.chdir(cwd)
    else:
        print(f"{pbs} not found, not submitting jobs for {family}.")

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
        log_message(message)
        raise argparse.ArgumentTypeError(message)
    else:
        if not os.path.getsize(filename) > 0:
            message = "{} file is empty".format(filename)
            print(message)
            raise argparse.ArgumentTypeError(message)
    return filename


if __name__ == "__main__":
    description = '''Reads sample info from TSV file (-f) and creates directory (-d) necessary to start
    crg2 from step (-s) requested.
    '''
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
            "-f",
            "--file",
            type=valid_file,
            required=True,
            help="Five column TAB-seperated sample info file",
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
    projects = {}
    Units = namedtuple("Units", "sample platform fq1 fq2 bam")
    platform = "ILLUMINA"
    args = parser.parse_args()
    with open(args.file) as f:
        for i in f:
            if not i.startswith("familyID"):
                family, sample, fq1, fq2, bam = i.strip("\n").split("\t")
                if not family in projects:
                    projects[family] = []
                projects[family].append(Units(sample, platform, fq1, fq2, bam))

    for i in projects:
        
        sample_list = projects[i]
        filepath = os.path.join(args.dir,i)

        #create project directory
        #copy config.yaml, dnaseq_cluster.pbs & replace family id; copy pbs_config.yaml 
        print(f"\nStarting to setup project directories for family: {i}")
        submit_flag = setup_directories(i, sample_list, args.dir, args.step)
        
        if submit_flag:
            write_proj_files(sample_list, filepath)
            submit_jobs(filepath, i)

    print("DONE")
