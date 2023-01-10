import argparse
import ast
import glob
import os
import sys
import pandas as pd
import pathlib
import re
import subprocess
from collections import namedtuple
from itertools import chain


"""
Usage: python3 genome_reanalysis.py -f <input_sample.tsv> -d <absolute path to create directories> 
Parses analyses request csv from Stager for genome re-analyses and sets up necessary directories (under 2nd argument), files as below:
1. create family analysis directory under directory provided
3. copy config_hpf.yaml, slurm-config.yaml and dnaseq_slurm_hpf.sh from crg2 repo and replace necessary strings
4. search results directory /hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes and /hpf/largeprojects/ccm_dccforge/dccdipg/c4r_wgs/results for crams from previous analyses
5. create units.tsv and samples.tsv for snakemake
6. submit job if all the above goes well
"""

crg2_dir = "/home/mcouse/crg2/"


def find_input(family, participant):
    """Given family and participant codenames, check if individuals have already been analyzed"""
    RESULTS_DIR_dccdipg = glob.glob(
        (f"/hpf/largeprojects/ccm_dccforge/dccdipg/c4r_wgs/results/{family}")
    )
    RESULTS_DIR_ccmmarvin = glob.glob(
        (f"/hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/{family}")
    )
    RESULTS_DIR = RESULTS_DIR_dccdipg + RESULTS_DIR_ccmmarvin
    try:
        RESULTS_DIR = pathlib.Path(RESULTS_DIR[0])
    except IndexError:
        RESULTS_DIR = None
    if RESULTS_DIR:
        print(RESULTS_DIR)
        cram_crg2 = RESULTS_DIR.rglob(f"{family}_{participant}.cram")
        cram_bcbio = RESULTS_DIR.rglob(f"{family}_{participant}-ready.cram")
        cram = chain(cram_crg2, cram_bcbio)
        print(f"{family}_{participant}.cram")
        try:
            input = str([c for c in cram][0])
        except IndexError:
            # family directory exists but cram does not
            input = None
    else:
        # family directory does not exist
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


def setup_directories(family, sample_list, filepath, analysis_id):
    # leave out analysis id for genomes for now
    d = os.path.join(filepath, family)
    if not os.path.isdir(d):
        os.makedirs(d)
        # copy config.yaml, pbs_config.yaml, and dnaseq_cluster.pbs
    for i in [
        "config_hpf.yaml",
        "slurm_profile/slurm-config.yaml",
        "dnaseq_slurm_hpf.sh",
    ]:
        cmd = ["cp", os.path.join(crg2_dir, i), d]
        subprocess.check_call(cmd)

    # replace family ID in config_hpf.yaml & dnaseq_slurm_cheo_ri.sh
    replace = f"s/NA12878/{family}/"
    pipeline = "s/wes/wgs/"
    PT_credentials = 's+PT_credentials: ""+PT_credentials: {}+'.format(
        "~/crg2/credentials.csv"
    )
    config = os.path.join(d, "config_hpf.yaml")
    if os.path.isfile(config):
        cmd = ["sed", "-i", replace, config]
        subprocess.check_call(cmd)
        cmd = ["sed", "-i", pipeline, config]
        subprocess.check_call(cmd)
        cmd = ["sed", "-i", PT_credentials, config]
        subprocess.check_call(cmd)
    replace = "s/job-name=crg2/job-name={}/".format(family)
    slurm = os.path.join(d, "dnaseq_slurm_hpf.sh")
    if os.path.isfile(slurm):
        cmd = ["sed", "-i", replace, slurm]
        subprocess.check_call(cmd)

    # glob hpo
    hpo_path = os.path.expanduser("~/gene_data/HPO")
    hpo = glob.glob(("{}/{}_HPO.txt").format(hpo_path, family))
    if len(hpo) > 1:
        print(f"Multiple HPO files found: {hpo}. Exiting!")
        exit()
    if len(hpo) == 1 and os.path.isfile(config):
        hpo = hpo[0]
        replace = 's+hpo: ""+hpo: "{}"+'.format(hpo)
        cmd = ["sed", "-i", replace, config]
        subprocess.check_call(cmd)

    #  write samples
    write_sample(filepath, family)

    # glob ped
    ped_path = os.path.expanduser("~/gene_data/pedigrees")
    ped = glob.glob(("{}/{}*ped").format(ped_path, family))
    if len(ped) > 1:
        print(f"Multiple ped files found: {ped}. Exiting!")
        exit()
    if len(ped) == 0:
        print(f"No ped files found for family: {family}")
        return False
    if len(ped) == 1 and os.path.isfile(config):
        ped = ped[0]
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
                samples = os.path.join(filepath, family, "samples.tsv")
                samples = pd.read_csv(samples, dtype=str)
                samples = list(samples.iloc[:, 0])
                samples = [family + s for s in samples]
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


def write_sample(filepath, family):
    samples = os.path.join(filepath, family, "samples.tsv")
    with open(samples, "w") as s:
        s.writelines("sample\n")
        for i in sample_list:
            s.writelines(f"{i.sample}\n")


def write_units(sample_list, filepath):
    units = os.path.join(filepath, "units.tsv")
    with open(units, "w") as u:
        u.writelines("sample\tplatform\tfq1\tfq2\tbam\tcram\n")
        for i in sample_list:
            u.writelines(
                f"{i.sample}\t{i.platform}\t{i.fq1}\t{i.fq2}\t{i.bam}\t{i.cram}\n"
            )


def submit_jobs(directory, family):
    if not os.path.isdir(directory):
        print(f"{directory} not found. Exiting!")
        exit()
    slurm = "sbatch"
    script = "dnaseq_slurm_cheo_ri.sh"
    cwd = os.getcwd()
    os.chdir(directory)
    cmd = [slurm, script]
    jobid = subprocess.check_output(cmd).decode("UTF-8").rstrip()
    print(f"Submitted snakemake job for {family}: {jobid}")
    os.chdir(cwd)


if __name__ == "__main__":
    description = "Reads sample info from Stager analysis csv file (-f) and creates directory (-d) necessary to run crg2"
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
        help="Absolute path where crg2 directory structure will be created under familyID/analysisID as base directory",
    )

    args = parser.parse_args()
    requested = pd.read_csv(args.file)
    inputs_list = []

    Units = namedtuple("Units", "sample platform fq1 fq2 bam cram")
    platform = "ILLUMINA"
    projects = {}
    analysis_ids = {}
    for index, row in requested.iterrows():
        family = list(set(ast.literal_eval(row["family_codenames"])))
        analysis_id = row["analysis_id"]
        # Will need to modify for multi-family analyses (i.e. cohort)
        family = family[0]
        if not family in projects:
            projects[family] = []
        analysis_ids[family] = analysis_id
        family_inputs = []
        for participant in ast.literal_eval(row["participant_codenames"]):
            # look for input files in results directory
            input = find_input(family, participant)
            print(input)
            if input:
                fq1, fq2, bam = "", "", ""
                cram = input
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
        # copy config_hpf.yaml, dnaseq_slurm_cheo_ri.sh & replace family id; copy slurm-config.yaml
        print(f"\nStarting to setup project directories for family: {i}")
        submit_flag = setup_directories(i, sample_list, args.dir, str(analysis_ids[i]))
        write_units(sample_list, filepath)

        if submit_flag:
            print(f"\nSubmitting job for family: {i}")
            submit_jobs(filepath, i)
        else:
            print(f"Job for {i} not submitted")
