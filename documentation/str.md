# ExpansionHunterDenovo version: 
    - location of the executable v0.7.0 is hard-coded in `config.yaml` file at
    config['tools']['ehdn']: "/hpf/largeprojects/ccmbio/arun/Tools/EHDN.TCAG/ExpansionHunterDenovo-v0.7.0". This version is not availble from GitHub or conda anymore, and the 1000G profiles from Brett Trost are based on this, so keeping this version as is.

# Scripts:
## R scripts:

1. scripts are fetched from `crg/str` repo
2. ~/crg/crg.ehdn.sh: runs the EHDN with many checks/searches for the relevant bam files based on run type: exome/genome and cre/crg/crg2 pipeline directory structures
3. ~/crg/ehdn_report.sh: main script that uses all the R packages from ~/crg/str, 1000G EHDN profiles, and generates script
    - ~/crg/str/DBSCAN.EHdn.parallel.R: library(dbscan), library(ggplot2), library(data.table), library(parallel), library(doSNOW)
    - ~/crg/str//mergeExpansions.R: library(data.table), library(GenomicRanges), library(ggplot2), library(cowplot), library(Biostrings)
    


## Python scripts:
1. ~/crg/str/compare_anchored_irrs.py: depends on "core" package local present
2. ~/crg/str/find_outliers.py: numpy, matplotlib
3. ~/crg/str/generate_EH_genotype_table.generic.py: argparse, glob, BTlib(local), collections, pandas, path
4. ~/crg/str/combine_counts.py: argparse, json, core(local)
5. ~/crg/str/BTlib.py: re, statistics, string, random, itertools, copy, collections, numpy, pandas, copy, docx
6. ~/crg/str/add_gene+threshold_to_EH_column_headings2.py: argparse
7. ~/crg/str/format_for_annovar.py: collections, re
8. ~/crg/str/format_from_annovar.py: pandas, xlsxwriter
    


