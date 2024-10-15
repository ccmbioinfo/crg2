#!/usr/bin/env python3

import argparse

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("EH_filename", type=str)
parser.add_argument("disease_locus_filename", type=str)
args = parser.parse_args()
#####################################

disease_locus_file = open(args.disease_locus_filename)

disease_locus_file.readline()

m = {}
for line in disease_locus_file:
    fields = line.rstrip("\n").split("\t")
    gene = fields[3]
    disease_threshold = fields[4]
    hg19_coords = fields[6].replace("chr", "")
    hg38_coords = fields[7]

    to_add = "{}:{}".format(gene, disease_threshold)
    if hg19_coords:
        m[hg19_coords] = to_add
    if hg38_coords:
        m[hg38_coords] = to_add


EH_file = open(args.EH_filename)
first_line = EH_file.readline()
fields = first_line.rstrip("\n").split("\t")
for i in range(0, len(fields)):
    for coords in m:
        if coords in fields[i]:
            if " (smaller allele)" in fields[i]:
                fields[i] = "{}:{}  (smaller allele)".format(fields[i].replace(" (smaller allele)", ""), m[coords])
            elif " (larger allele)" in fields[i]:
                    fields[i] = "{}:{}  (larger allele)".format(fields[i].replace(" (larger allele)", ""), m[coords])
            else:
                fields[i] = "{}:{}".format(fields[i], m[coords])
print("\t".join(fields))
for line in EH_file:
    print(line, end="")
