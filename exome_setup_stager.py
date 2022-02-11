import argparse
import ast
import logging
import os
import shutil
import sys

# replace with /home/ccmmarvin/crg2 when in production
crg2_dir = os.path.dirname(os.path.realpath(sys.argv[0]))


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
    for i in ["config.yaml", "pbs_profile/pbs_config.yaml", "dnaseq_cluster.pbs"]:
        src = os.path.join(crg2_dir, i)
        if i == "pbs_profile/pbs_config.yaml":
            dest = os.path.join(filepath, "pbs_config.yaml")
        else:
            dest = os.path.join(filepath, i)
        shutil.copyfile(src, dest)


def write_config(family, input_type_list, filepath):
    # replace family ID in config.yaml & dnaseq_cluster.pbs
    config = os.path.join(filepath, "config.yaml")
    replace_str(config, "NA12878", family)

    # modify input type depending on input file
    input_set = set(input_type_list)
    if input_set == {"bam"}:
        replace_str(config, 'input: "fastq"', 'input: "bam"')
    elif input_set == {"cram"}:
        replace_str(config, 'input: "fastq"', 'input: "cram"')
    elif input_set != {"fastq"}:
        logging.error("Mixed input filetypes among participants detected, exiting")
        sys.exit(1)


def input_type(file, participant):
    # determine filetype of linked file based on file suffix
    if file.endswith("fastq.gz") or file.endswith("fq.gz"):
        if "R1" in file:
            filetype = "fq1"
        elif "R2" in file:
            filetype = "fq2"
        else:
            logging.error(
                "Unrecognized fastq file format for %s, exiting" % participant
            )
            sys.exit(1)
    elif file.endswith("bam"):
        filetype = "bam"
    elif file.endswith("cram"):
        filetype = "cram"
    else:
        filetype = None
    return filetype


def dataset_to_dict(datasets):
    # create dict of dicts with participant/sample name and associated files
    all_datasets_dict = {}
    for dataset in datasets:
        participant = ast.literal_eval(dataset.split(":")[0])
        files = ast.literal_eval(dataset.split(":")[1])
        keys = ["fq1", "fq2", "bam", "cram"]
        dataset_dict = {k: [] for k in keys}
        all_datasets_dict[participant] = dataset_dict
        for f in files:
            filetype = input_type(f, participant)
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
    input_type_list = []
    with open(units, "w") as u, open(samples, "w") as s:
        s.writelines("sample\n")
        u.writelines("sample\tplatform\tfq1\tfq2\tbam\tcram\n")
        for participant in datasets_dict:
            s.writelines(f"{participant}\n")
            bam = datasets_dict[participant]["bam"]
            cram = datasets_dict[participant]["cram"]
            fq1 = datasets_dict[participant]["fq1"]
            fq2 = datasets_dict[participant]["fq2"]
            # if fastqs available, set these as input files
            if len(fq1) != 0:
                fq1 = ",".join(fq1)
                fq2 = ",".join(fq2)
                u.writelines(f"{participant}\tILLUMINA\t{fq1}\t{fq2}\t\t\n")
                input_type_list.append("fastq")
            # if bam available, set as input file
            elif len(bam) != 0:
                if len(bam) == 1:
                    u.writelines(f"{participant}\tILLUMINA\t\t\t{bam[0]}\t\n")
                    input_type_list.append("bam")
                else:
                    logging.error(
                        "Multiple bam files provided for %s, exiting" % participant
                    )
                    sys.exit(1)
            # if cram available, set as input file
            elif len(cram) != 0:
                if len(cram) == 1:
                    u.writelines(f"{participant}\tILLUMINA\t\t\t{cram[0]}\t\n")
                    input_type_list.append("cram")
                else:
                    logging.error("Multiple cram files provided, exiting")
                    sys.exit(1)
            else:
                logging.error("No input file provided for %s, exiting" % participant)
                sys.exit(1)
    return input_type_list


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
        nargs="+",
        required=True,
        help="List of participants and their associated files in the format participant_codename:linked_file(s)",
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
    filepath = os.path.join(args.analysis_directory, family)

    logging.basicConfig(
        filename="%s/setup.log" % filepath,
        level=logging.INFO,
        format="%(asctime)s %(message)s",
    )

    logging.info("Setting up directory for crg2 run")
    setup_directory(filepath)
    logging.info("Parsing participants' linked files")
    datasets_dict = dataset_to_dict(datasets)
    logging.info("Writing units.tsv and samples.tsv files")
    input_type_list = write_units_samples(datasets_dict, filepath)
    logging.info("Parsing config.yaml")
    write_config(family, input_type_list, filepath)
    logging.info("crg2 set up complete")
