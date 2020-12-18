import sys, os, subprocess
from collections import namedtuple

sample_sheet = sys.argv[1]
path = sys.argv[2]
crg2_dir = os.path.realpath(sys.argv[0])

def create_symlink(src, dest):
    if not os.path.isfile(dest):
        print(f"symlinking {src} -> {dest}")
        #subprocess.checkcall(["ln", "-s", src, dest])

def setup_directories(s, sample_list, filepath)
    d = os.path.join(filepath, family)
    if not os.path.isdir(d):
        print(f"making dir {d}")
        #os.mkdir(d)
    #copy config.yaml, pbs_config.yaml, and dnaseq_cluster.pbs
    for i in ["config.yaml", "pbs_profile/pbs_config.yaml", "dnaseq_cluster.pbs"]
        cmd = ["cp", os.path.join(crg2_dir, i), d]
        print(f"copying {cmd}")
        # subprocess.check_call(cmd)
        
    #bam: start after remove_decoy 
    if len(sample_list) == len([i for i in sample_list if i.bam ]):
        decoy = os.path.join(d, "decoy_rm")
        print(f"creating dir {deco}")
        # os.mkdir(decoy)
        for i in sample_list:
            dest = "{}-{}.no_decoy_reads.bam".format(i.sample,i.unit)
            create_symlink(i.bam, os.path.join(decoy, dest))
            if os.path.isfile(i.bam + ".bai"):
                dest = dest + ".bai"
                create_symlink(src, os.path.join(decoy, dest))

def write_proj_files(sample_list, filepath):
    units, samples = os.path.join(filepath, "units.tsv"), os.path.join(filepath, "samples.tsv")
    with open(units, "w") as u and open(samples) as s:
        u.writelines("sample\n")
        s.writelines("sample\tunit\tplatform\tfq1\tfq2\tbam\n")
        for i in sample_list:
            s.writelines(f"{i.sample}\t{i.unit}\t{i.platform}\{i.fq1}\t{i.fq2}\t{i.bam}\n")
            u.writelines(f"{i.sample}\n")                

projects = {}
Units = namedtuple("Units", "sample unit platform fq1 fq2 bam")
unit, platform = 1, "ILLUMINA"
with open(sample_sheet) as f:
    for i in f:
        family, sample, inp = i.strip().split("\t")
        bam, fq1, fq1 = "", "", ""
        if ".bam" in inp:
            bam = inp
        elif ".fastq" in inp:
            if "," in inp:
                #multiple fastq
                pass
            else:
                bpath = os.path.dirname(inp)
                prefix, suffix = os.path.basename(inp).split("_R1")
                r2 = os.path.join(bpath, prefix + "_R2" + suffix)
                if os.path.isfile(r2) and os.path.isfile(inp):
                    fq1, fq2 = inp, r2
                else:
                    print(f"{inp} or {r2} was not found")
        else:
            print("unhandled input type! exiting")
            exit()
        if not family in projects:
            projects[family] = []
        projects[family].append(Units(sample, unit, platform, fq1, fq2, bam))

for i in projects:
    
    sample_list = projects[i]
    filepath = os.path.join(path,i)

    #create project directory
    setup_directories(i, sample_list, path)
    
    #create samples.tsv & units.tsv
    write_proj_files(sample_list, filepath)
    
    #copy config.yaml & replace project: 
    #copy dnaseq_cluster.pbs & pbs_config.yaml 






