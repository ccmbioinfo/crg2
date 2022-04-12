from posixpath import realpath
import sys, os
from pathlib import Path


def read_idmap(idmap_file):
    sample_map = {}
    with open(idmap_file) as f:
        for i in f:
            sample, id = i.strip().split("\t")
            sample_map[id] = sample
    return sample_map


if len(sys.argv) <= 1:
    print(
        "Usage: python gen_sample_sheet.py </path/to/download_folder> [optional:map_id file]"
    )
    print(
        "Output: 2 TSV files named <folder_name>-samples_for_crg2-new.tsv and <folder_name>-samples_for_crg.tsv"
    )
    print("Required argument: folder name was not passed. Exiting!")
    print(
        "If FASTQ files are not named by sample names, please pass the 2nd argument with file id -> sample mapping"
    )
    sys.exit()

folder_name = sys.argv[1]
projects = {}

if len(sys.argv) > 2:
    idmap_file = sys.argv[2]
    idmap = read_idmap(idmap_file)
    projects = {}
    for fastq_id in idmap:
        sampleid = idmap[fastq_id]
        if not sampleid in projects:
            projects[sampleid] = {"read1": "", "read2": ""}
        projects[sampleid]["read1"] = [
            str(os.path.realpath(path))
            for path in Path(folder_name).rglob(fastq_id + "*R1*f*q.gz")
        ]
        print(
            fastq_id,
            [
                str(os.path.realpath(path))
                for path in Path(folder_name).rglob(fastq_id + "*R1*f*q.gz")
            ],
        )
        # print(projects[sampleid]["read1"], fastq_id, sampleid)
        projects[sampleid]["read2"] = [
            str(os.path.realpath(path))
            for path in Path(folder_name).rglob(fastq_id + "*R2*f*q.gz")
        ]

else:
    read1 = [path for path in Path(folder_name).rglob("*R1*fastq.gz")]
    read2 = [path for path in Path(folder_name).rglob("*R2*fastq.gz")]
    for i in read1:
        key = os.path.basename(i).split("_RGS_")[0]
        print(key)
        if not key in projects:
            projects[key] = {"read1": [], "read2": []}
        projects[key]["read1"].append(str(os.path.realpath(i)))

    for i in read2:
        key = os.path.basename(i).split("_RGS_")[0]
        print(key)
        if not key in projects:
            projects[key] = {"read1": [], "read2": []}
        projects[key]["read2"].append(str(os.path.realpath(i)))


with open(folder_name + "-samples_for_crg2-new.tsv", "w") as f, open(
    folder_name + "-samples_for_crg.tsv", "w"
) as f1:
    f1.writelines(f"sampleid\tfastq\tbam\n")
    for i in sorted(projects):
        project, sample = i.split("_")
        # sort both ends to match read order while concatenating; else would fail during mapping
        fq1 = ",".join(sorted(projects[i]["read1"]))
        fq2 = ",".join(sorted(projects[i]["read2"]))
        # for crg2
        f.writelines(f"{project}\t{sample}\t{fq1}\t{fq2}\t\n")
        # for crg
        f1.writelines(f"{i}\t{fq1}\t\n")
