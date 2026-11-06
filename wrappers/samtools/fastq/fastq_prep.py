import pandas as pd
import sys, os, subprocess, logging, re


def log_message(*message):
    """write message to logfile and stdout"""
    if message:
        for i in message:
            logging.info(i)
            print(i)


def create_symlink(src, dest):
    """symlink file"""
    if not os.path.isfile(dest):
        log_message(f"creating symlinks {src} to {dest}")
        subprocess.check_call(["ln", "-s", src, dest])


def find_read_id(filename, read_number):
    fq_regex = re.compile(f"_{read_number}")
    fastqname, dirname = os.path.basename(filename), os.path.dirname(filename)
    rid = re.findall(fq_regex, fastqname)
    if rid:
        split_str = rid[0]
    else:
        log_message(
            'Unhandled fastq read identifier. Read identifier should be "_R1" or "_R2". Exiting! '
        )
        exit()
    return fastqname, dirname, split_str


def check_fastq(read1, read2):  # multiple fastq files per end
    """make sure that each R1 has a corresponding R2 and vice versa"""
    for read in read1:
        fastqname, dirname, split_str = find_read_id(read, "R1")
        prefix, suffix = fastqname.split(split_str)
        paired_read = os.path.join(dirname, prefix + "_R2" + suffix)
        if paired_read not in read2:
            log_message(f"{paired_read} not in {read2}")
            exit()
    for read in read2:
        fastqname, dirname, split_str = find_read_id(read, "R2")
        prefix, suffix = fastqname.split(split_str)
        paired_read = os.path.join(dirname, prefix + "_R1" + suffix)
        if paired_read not in read1:
            log_message(f"{paired_read} not in {read1}")
            exit()


def concatenate_fastq(r1, r2, family, sample):
    """
    based on input file suffix, performs
    1. concatenation (more than 1 fastq per end)
    2. create symlink if only one fastq per end
    """
    if len(r1) > 1:  # multiple fastq files per end
        log_message(f"Multiple fastq files per end {sample}")
        # check_fastq(r1, r2)
        r1_args = " ".join(r1)
        r2_args = " ".join(r2)
        cmd_r1 = [f"cat {r1_args} > fastq/{family}_{sample}_R1.fastq.gz"]
        cmd_r2 = [f"cat {r2_args} > fastq/{family}_{sample}_R2.fastq.gz"]
        log_message(f"Command: {cmd_r1}")
        subprocess.run(cmd_r1, shell=True)
        log_message(f"Command: {cmd_r2}")
        subprocess.run(cmd_r2, shell=True)

    elif ".fastq" in r1[0] or ".fq" in r1[0] and ".fastq" in r2[0] or ".fq" in r2[0]:
        # one fastq per end; symlink
        log_message(f"Single FASTQ file per end for {sample} ")
        # check_fastq(r1, r2)
        dest = os.path.join(f"fastq/{family}_{sample}_R1.fastq.gz")
        create_symlink(r1[0], dest)
        dest = os.path.join(f"fastq/{family}_{sample}_R2.fastq.gz")
        create_symlink(r2[0], dest)
    else:
        print(r1, r2)
        log_message(f"Input {r1}, {r2} given for {sample} is not handled. Exiting!")
        exit()

def decompress_ora(r1,r2,orad,orad_ref):
    decompressed_r1_list=[]
    decompressed_r2_list=[]
    for read in r1:
        log_message(f"{read}")
        command = f"{orad} --ora-reference {orad_ref} {read}"
        log_message(f"Command:{command}")
        run_orad=subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        log_message(run_orad.stdout)
        #log_message(run_orad.stderr)
        decompressed_r1_list.append(f"{os.path.splitext(read)[0]}.gz")
    
    log_message(f"R1 List: {decompressed_r1_list}")

    for read in r2:
        log_message(f"{read}")
        command = f"{orad} --ora-reference {orad_ref} {read}"
        log_message(f"Command:{command}")
        run_orad=subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        log_message(run_orad.stdout)
        #log_message(run_orad.stderr)
        decompressed_r2_list.append(f"{os.path.splitext(read)[0]}.gz")

    log_message(f"R2 List: {decompressed_r2_list}")
    
    return decompressed_r1_list,decompressed_r2_list

def remove_bgzipped_files(decompressed_r1_list,decompressed_r2_list):
    if len(decompressed_r1_list+decompressed_r2_list) > 2:
        for read in decompressed_r1_list+decompressed_r2_list:
            log_message(f"Removing {read}")
            command=f"rm {read}"
            subprocess.run(command,shell=True)
    else:
        log_message(f"No files to remove")

def main(units, sample, family, orad, orad_ref):
    logfile = f"logs/input_prep/{family}_{sample}.log"
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.DEBUG,
        format="%(asctime)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M",
    )
    # extract read1, read2 from units.tsv
    fastqs = units.loc[(sample), ["fq1", "fq2"]].dropna()
    r1, r2 = fastqs.fq1, fastqs.fq2
    r1, r2 = r1.split(","), r2.split(",")
    # check if fastq files exist
    for read in r1 + r2:
        if not os.path.isfile(read):
            log_message(f"{read} does not exist. Exiting!")
            exit()
    # concatenate fastqs, or symlink if only one per end

    #Checking if the fastq file extension is .gz or .ora
    if os.path.splitext(r1[0])[1] == ".gz":
        #If file is gzipped then concatenate the reads
        log_message(f"Files for {sample} are bgzipped.")
        concatenate_fastq(r1, r2, family, sample)
    elif os.path.splitext(r1[0])[1] == ".ora":
        #If file is ora compressed, then first decompress it
        log_message(f"Files for {sample} are in ora format. Decompressing the files using orad")
        decompressed_r1, decompressed_r2=decompress_ora(r1,r2,orad,orad_ref)
        concatenate_fastq(decompressed_r1, decompressed_r2, family, sample)
        #remove the bgzipped files
        remove_bgzipped_files(decompressed_r1,decompressed_r2)


if __name__ == "__main__":
    units = snakemake.input.units
    sample = snakemake.wildcards.sample
    family = snakemake.wildcards.family
    orad = snakemake.wildcards.orad
    orad_ref = snakemake.wildcards.orad_ref
    main(units, sample, family, orad, orad_ref)
