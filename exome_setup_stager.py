import argparse
import ast
from datetime import datetime
import json
import logging
import os
import shutil
import sys
import re

crg2_dir = "/srv/shared/pipelines/crg2"


def replace_str(filename, target, replacement):
    # pythonized sed
    with open(filename, "r") as file:
        lines = file.readlines()
    with open(filename, "w") as file:
        for line in lines:
            line = line.replace(target, replacement)
            file.write(line)


def setup_directory(filepath):
    # copy config.yaml, pbs_config.yaml, and dnaseq_cluster.pbs
    for i in ["config.yaml", "slurm_profile/slurm-config.yaml", "dnaseq_slurm.sh"]:
        src = os.path.join(crg2_dir, i)
        if i == "slurm_profile/slurm-config.yaml":
            dest = os.path.join(filepath, "slurm-config.yaml")
        else:
            dest = os.path.join(filepath, i)
        shutil.copyfile(src, dest)


def write_config(family, filepath):
    # replace family ID in config.yaml & dnaseq_cluster.pbs
    config = os.path.join(filepath, "config.yaml")
    replace_str(config, "NA12878", family)


def input_type(file, participant):
    # determine filetype of linked file based on file suffix
    if file.endswith("fastq.gz") or file.endswith("fq.gz"):
        filetype = fastq_regex(file)
    elif file.endswith("bam"):
        filetype = "bam"
    elif file.endswith("cram"):
        filetype = "cram"
    else:
        filetype = None
    return filetype

def fastq_regex(fastq):
    # regex would capture reads such as 19-13210-A-02-00_BH2HKGBCX3_R1_1.fastq.gz, or 670513-P_HVWLCBCXX-2-ID01_1_sequence.fastq.gz
    read1_regex = re.compile(r"_R1[_.]|_1_")
    read2_regex = re.compile(r"_R2[_.]|_2_")
    match1 = re.search(read1_regex, fastq)
    match2 = re.search(read2_regex, fastq)
    if match1:
        return "fq1"
    elif match2:
        return "fq2"
    else:
        logging.error(
            "Unrecognized fastq file format for %s, exiting" % participant
        )
        sys.exit(1)

def dataset_to_dict(datasets):
    # create dict of dicts with participant/sample name and associated files
    datasets = json.loads(datasets)
    all_datasets_dict = {}
    for participant in datasets:
        files = datasets[participant]
        logging.info("creating dataset dictionary for participant %s", participant)
        keys = ["fq1", "fq2", "bam", "cram"]
        dataset_dict = {k: [] for k in keys}
        all_datasets_dict[participant] = dataset_dict
        for f in files:
            filetype = input_type(f, participant)
            f = os.path.join('/srv/minio', f)
            logging.info("adding file %s",  f)
            if filetype == "bam":
                all_datasets_dict[participant]["bam"].append(f)
            elif filetype == "cram":
                all_datasets_dict[participant]["cram"].append(f)
            elif filetype == "fq1":
                all_datasets_dict[participant]["fq1"].append(f)
            elif filetype == "fq2":
                all_datasets_dict[participant]["fq2"].append(f)
    return all_datasets_dict


def write_units_samples(datasets_dict, filepath):
    units, samples = os.path.join(filepath, "units.tsv"), os.path.join(
        filepath, "samples.tsv"
    )
    with open(units, "w") as u, open(samples, "w") as s:
        s.writelines("sample\n")
        u.writelines("sample\tplatform\tfq1\tfq2\tbam\tcram\n")
        for participant in datasets_dict:
            s.writelines(f"{participant}\n")
            bam = datasets_dict[participant]["bam"]
            cram = datasets_dict[participant]["cram"]
            fq1 = sorted(datasets_dict[participant]["fq1"])
            print(fq1)
            print(datasets_dict)
            fq2 = sorted(datasets_dict[participant]["fq2"])
            # if fastqs available, set these as input files
            if len(fq1) != 0:
                fq1 = ",".join(fq1)
                fq2 = ",".join(fq2)
                u.writelines(f"{participant}\tILLUMINA\t{fq1}\t{fq2}\t\t\n")
            # if bam available, set as input file
            elif len(bam) != 0:
                if len(bam) == 1:
                    u.writelines(f"{participant}\tILLUMINA\t\t\t{bam[0]}\t\n")
                else:
                    logging.error(
                        "Multiple bam files provided for %s, exiting" % participant
                    )
                    sys.exit(1)
            # if cram available, set as input file
            elif len(cram) != 0:
                if len(cram) == 1:
                    u.writelines(f"{participant}\tILLUMINA\t\t\t\t{cram[0]}\n")
                else:
                    logging.error("Multiple cram files provided, exiting")
                    sys.exit(1)
            else:
                logging.error("No input file provided for %s, exiting" % participant)
                sys.exit(1)



if __name__ == "__main__":
    description = "Sets up crg2 exome pipeline run from Stager analysis metadata"
    parser = argparse.ArgumentParser(description=description)
    # set up directory:
    # copy config files etc
    # based on family, pipeline and input type, modify config.yaml
    # create units.tsv and samples.tsv based on input
    parser.add_argument(
        "-f",
        "--family",
        type=str,
        required=True,
        help="Family codename",
    )

    parser.add_argument(
        "-d",
        "--datasets",
        type=str,
        required=True,
        help="json string listing participants and their associated files",
    )

    parser.add_argument(
        "-a",
        "--analysis_directory",
        type=str,
        required=True,
        help="Base path for in which analyses are run",
    )

    args = parser.parse_args()
    family = args.family
    datasets = args.datasets
    filepath = args.analysis_directory
    today = datetime.today().strftime('%Y-%m-%d')

    logging.basicConfig(
        filename="%s/setup_%s.log" % (filepath, today),
        level=logging.INFO,
        format="%(asctime)s %(message)s",
    )

    logging.info("Setting up directory for crg2 run")
    setup_directory(filepath)
    logging.info("Parsing participants' linked files")
    datasets_dict = dataset_to_dict(datasets)
    logging.info("Writing units.tsv and samples.tsv files")
    write_units_samples(datasets_dict, filepath)
    logging.info("Parsing config.yaml")
    write_config(family, filepath)
    logging.info("crg2 set up complete")
