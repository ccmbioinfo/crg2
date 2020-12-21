import sys, os, subprocess
from collections import namedtuple

'''
Usage: python parser.py <input_sample.tsv> <absolute path to create directories> 
Parse three-column(family,sample,file-path) TSV file (1st argument) and sets up necessary directories (under 2nd argument), files as below:
1. create family and decoy_rm directory
2. symlink BAM files to deocy_rm
3. copy config.yaml, pbs_config.yaml and dnaseq_cluster.pbs from crg2 repo and replace necessary string
4. create units.tsv and samples.tsv for snakemake
5. submit job if all the above goes well
'''

sample_sheet = sys.argv[1]
path = sys.argv[2]
crg2_dir = os.path.dirname(os.path.realpath(sys.argv[0]))

def create_symlink(src, dest):
    if not os.path.isfile(dest):
        cmd = ["ln", "-s", src, dest]
        subprocess.check_call(cmd)

def setup_directories(family, sample_list, filepath):
    d = os.path.join(filepath, family)
    if not os.path.isdir(d):
        os.mkdir(d)

    #copy config.yaml, pbs_config.yaml, and dnaseq_cluster.pbs
    for i in ["config.yaml", "pbs_profile/pbs_config.yaml", "dnaseq_cluster.pbs"]:
        cmd = ["cp", os.path.join(crg2_dir, i), d]
        subprocess.check_call(cmd)

    #replace family ID in config.yaml & dnaseq_cluster.pbs
    replace = 's/NA12878/{}/'.format(family)
    config = os.path.join(d,"config.yaml")
    if os.path.isfile(config):
        cmd = ["sed", "-i", replace, config]
        subprocess.check_call(cmd)
    replace = 's/crg2_pbs/{}/'.format(family)
    pbs = os.path.join(d,"dnaseq_cluster.pbs")
    if os.path.isfile(pbs):
        cmd = ["sed", "-i", replace, pbs]
        subprocess.check_call(cmd)

    #bam: start after remove_decoy 
    if len(sample_list) == len([i for i in sample_list if i.bam ]):
        decoy = os.path.join(d, "decoy_rm")
        if not os.path.isdir(decoy):
            os.mkdir(decoy)
            for i in sample_list:
                dest = "{}-{}.no_decoy_reads.bam".format(i.sample,i.unit)
                create_symlink(i.bam, os.path.join(decoy, dest))
                if os.path.isfile(i.bam + ".bai"):
                    src = i.bam + ".bai"
                    dest = dest + ".bai"
                    create_symlink(src, os.path.join(decoy, dest))
        else:
            print(f"{decoy} directory already exists! Won't rewrite/submit jobs. ")
            return False
    return True


def write_proj_files(sample_list, filepath):
    units, samples = os.path.join(filepath, "units.tsv"), os.path.join(filepath, "samples.tsv")
    with open(units, "w") as u, open(samples, "w") as s:
        s.writelines("sample\n")
        u.writelines("sample\tunit\tplatform\tfq1\tfq2\tbam\n")
        for i in sample_list:
            s.writelines(f"{i.sample}\n")
            u.writelines(f"{i.sample}\t{i.unit}\t{i.platform}\t{i.fq1}\t{i.fq2}\t{i.bam}\n")
            

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


projects = {}
Units = namedtuple("Units", "sample unit platform fq1 fq2 bam")
unit, platform = 1, "ILLUMINA"
with open(sample_sheet) as f:
    for i in f:
        if not i.startswith("familyID"):
            family, sample, inp = i.strip().split("\t")
            bam, fq1, fq2 = "", "", ""
            if ".bam" in inp:
                if os.path.isfile(inp):
                    bam = inp
                else:
                    print(f"\tBAM: {inp} was not found! Exiting!")
                    exit()
            elif ".fastq" in inp:
                if "," in inp:
                    #multiple fastq
                    pass
                else:
                    pass
                    # bpath = os.path.dirname(inp)
                    # prefix, suffix = os.path.basename(inp).split("_R1")
                    # r2 = os.path.join(bpath, prefix + "_R2" + suffix)
                    # if os.path.isfile(r2) and os.path.isfile(inp):
                    #     fq1, fq2 = inp, r2
                    # else:
                    #     print(f"{inp} or {r2} was not found")
            else:
                print("Unhandled input type! Exiting")
                exit()
            if not family in projects:
                projects[family] = []
            projects[family].append(Units(sample, unit, platform, fq1, fq2, bam))

for i in projects:
    
    sample_list = projects[i]
    filepath = os.path.join(path,i)

    #create project directory
    #copy config.yaml, dnaseq_cluster.pbs & replace family id; copy pbs_config.yaml 
    print(f"\nStarting to setup project directories for family: {i}")
    submit_flag = setup_directories(i, sample_list, path)
    
    if submit_flag:
        write_proj_files(sample_list, filepath)
        submit_jobs(filepath, i)

print("DONE")

