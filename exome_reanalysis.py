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


"""
Usage: python3 exome_reanalysis.py -f <input_sample.tsv> -d <absolute path to create directories> 
Parses analyses request csv from Stager for exome re-analyses and sets up necessary directories (under 2nd argument), files as below:
1. create family analysis directory under directory provided
3. copy config_cheo_ri.yaml, slurm-config.yaml and dnaseq_slurm_cheo_ri.sh from crg2 repo and replace necessary strings
4. search results directory /srv/shared/hpf/exomes/results for crams from previous analyses
5. create units.tsv and samples.tsv for snakemake
6. submit job if all the above goes well
"""

crg2_dir = "/srv/shared/pipelines/crg2/"
aliases=f"{crg2_dir}/fam_ptp-2022-07-11.tsv"

def create_alias_long(stager_df:pd.DataFrame, main_col = 'participant_codename',
                      alias_col = 'participant_aliases') -> pd.DataFrame:
    """From Delvin So's G4RD utils script"""
    ss_aliases = pd.concat([stager_df[main_col],
                        stager_df[alias_col].str.split(r',|;',expand=True, regex = True)], axis = 1)\
    .rename(lambda x: 'alias' + str(x) if len(str(x)) == 1 else x, axis = 1)


    aliases_long = pd.melt(ss_aliases.reset_index(),
        id_vars = [main_col],
        value_vars = [x for x in ss_aliases.columns if str(x).startswith('alias')])

    aliases_long['value'] = aliases_long['value'].str.strip().drop_duplicates()
    aliases_long = aliases_long.dropna(subset = ['value'])

    same_df = pd.DataFrame.from_dict({main_col : stager_df[main_col].drop_duplicates().tolist(),
                        'variable': ['codename' for x in range(len(stager_df[main_col].drop_duplicates().tolist()))],
                        'value' : stager_df[main_col].drop_duplicates().tolist()})

    aliases_long_final = pd.concat([aliases_long, same_df], axis = 0)
    aliases_long_final['type'] = main_col.split('_')[0]

    return aliases_long_final

def check_alias(participant, alias_df):
    try: 
        print(participant)
        mapping = alias_df[alias_df["participant_codename"] == participant]["value"].values[0]
    except:
        mapping = None
    return mapping 

def parse_aliases(aliases):
    stager_df = pd.read_csv(aliases, delimiter = '\t',encoding='cp1252')
    stager_df_ss = stager_df[['participant_id', 'family_id', 'participant_codename', 'participant_aliases', 'participant_type', 'family_codename', 'family_aliases']]
    ptp_aliases_long = create_alias_long(stager_df = stager_df_ss)
    return ptp_aliases_long

def find_input(family, participant):
    """Given family and participant codenames, check if individuals have already been analyzed"""
    RESULTS_DIR = glob.glob((f"/srv/shared/hpf/exomes/results/*/{family}*"))
    try:
        RESULTS_DIR = pathlib.Path(RESULTS_DIR[0])
    except IndexError:
        RESULTS_DIR = None
    if RESULTS_DIR:
        family = os.path.basename(RESULTS_DIR)
    if RESULTS_DIR:
        print(RESULTS_DIR)
        cram = RESULTS_DIR.rglob(f"{family}_{participant}.cram")
        try:
            input = str([c for c in cram][0])
        except IndexError:
            # family directory exists but cram does not
            # bam file may be named under an alias
            alias_mapping = parse_aliases(aliases)
            mapping = check_alias(participant, alias_mapping)
            cram = RESULTS_DIR.rglob(f"{family}_{mapping}.cram")
            cram = [c for c in cram]
            if len(cram) == 0:
                input = None
            else:
                input = str(cram[0])
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
    d = os.path.join(filepath, family, analysis_id)
    if not os.path.isdir(d):
        os.makedirs(d)
        # copy config.yaml, pbs_config.yaml, and dnaseq_cluster.pbs
    for i in [
        "config_cheo_ri.yaml",
        "slurm_profile/slurm-config.yaml",
        "dnaseq_slurm_cheo_ri.sh",
    ]:
        cmd = ["cp", os.path.join(crg2_dir, i), d]
        subprocess.check_call(cmd)

    # replace family ID in config_hpf.yaml & dnaseq_slurm_cheo_ri.sh
    replace = f"s/NA12878/{family}/"
    config = os.path.join(d, "config_cheo_ri.yaml")
    if os.path.isfile(config):
        cmd = ["sed", "-i", replace, config]
        subprocess.check_call(cmd)
    replace = f"s/crg2_pbs/{family}/"
    pbs = os.path.join(d, "dnaseq_slurm_cheo_ri.sh")
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
        filepath = os.path.join(args.dir, i, str(analysis_ids[i]))

        # create project directory
        # copy config_hpf.yaml, dnaseq_slurm_cheo_ri.sh & replace family id; copy slurm-config.yaml
        print(f"\nStarting to setup project directories for family: {i}")
        submit_flag = setup_directories(i, sample_list, args.dir, str(analysis_ids[i]))

        write_proj_files(sample_list, filepath)
        if submit_flag:
            print(f"\nSubmitting job for family: {i}")
            submit_jobs(filepath, i)
        else:
            print(f"Job for {i} not submitted")
