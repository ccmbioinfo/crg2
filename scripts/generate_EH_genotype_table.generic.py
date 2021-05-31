#!/usr/bin/env python3

import argparse
import BTlib
import glob
import os
from collections import defaultdict
import pandas
from pathlib import Path

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("base_dir", type=str)
parser.add_argument("--two-column-per-locus-format", default=False, action="store_true")
parser.add_argument("--skip-gender", default=False, action="store_true") # Skip VCF files containing "male" or "female" in the filename (specific to PPMI data)
args = parser.parse_args()
#####################################

mat = defaultdict(lambda: defaultdict(lambda: "./."))

all_samples = set()
all_locus_keys = {}
total_files_processed = 1

for vcf_filename in Path(args.base_dir).rglob('*.vcf'):
    if args.skip_gender and "male" in str(vcf_filename): # Skip VCF files containing "male" or "female" in the filename (specific to PPMI data)
        continue

    BTlib.eprint("Now processing file {} ({})".format(total_files_processed, vcf_filename))
    total_files_processed += 1
    vcf_file = open(vcf_filename)

    # Figure out name of sample
    sample = os.path.basename(vcf_filename).split(".")[0].replace("_recal-mini", "")
    all_samples.add(sample)

    for line in vcf_file:
        if line[0] == "#":
            continue
        fields = line.split("\t")
        chrom = fields[0]
        start = fields[1]
        end = fields[7].split(";")[0].split("=")[1]
        repeat_unit = fields[7].split(";")[3].split("=")[1]

        locus_key = "{}:{}-{}:{}".format(chrom, start, end, repeat_unit)
        all_locus_keys[locus_key] = True

        genotypes = [] # 0 = reference size, 1 = first alternate size, etc.
        genotypes.append(int(fields[7].split(";")[1].split("=")[1])) # Get reference size

        alt = fields[4].split(",") # Get list of alternate sizes

        for i in range(0, len(alt)):
            if alt[i] != ".":
                genotypes.append(int(alt[i].replace("<STR", "").replace(">", "")))

        gt = fields[9].split(":")[0]

        if "/" in gt: # Two alleles
            allele1, allele2 = [int(x) for x in gt.split("/")]
            min_size = min(genotypes[allele1], genotypes[allele2])
            max_size = max(genotypes[allele1], genotypes[allele2])
            mat[sample][locus_key] = "{}/{}".format(min_size, max_size)
        else: # One allele (X in males)
            mat[sample][locus_key] = str(genotypes[int(gt)])

all_samples_sorted = sorted(all_samples)
all_locus_keys_sorted = sorted(all_locus_keys)
if args.two_column_per_locus_format:
    print("Sample", end="")
    for locus_key in all_locus_keys_sorted:
        print("\t{} (smaller allele)\t{} (larger allele)".format(locus_key, locus_key), end="")
    print()
else:
    print("Sample\t{}".format("\t".join(all_locus_keys_sorted)))
for sample in all_samples_sorted:
    if args.two_column_per_locus_format:
        print(sample, end="")
        for locus_key in all_locus_keys_sorted:
            if "/" in mat[sample][locus_key]:
                smaller_allele = mat[sample][locus_key].split("/")[0]
                larger_allele = mat[sample][locus_key].split("/")[1]
            else:
                smaller_allele = "."
                larger_allele = mat[sample][locus_key]

            print("\t{}\t{}".format(smaller_allele, larger_allele), end="")
        print()
    else:
        print("{}\t{}".format(sample, "\t".join([mat[sample][locus_key] for locus_key in all_locus_keys_sorted])))
