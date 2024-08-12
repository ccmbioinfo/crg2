ExpansionHunterDenovo: location of the executable v0.7.0 is hard-coded in `config.yaml` file at
    config['tools']['ehdn']: "/hpf/largeprojects/ccmbio/arun/Tools/EHDN.TCAG/ExpansionHunterDenovo-v0.7.0"
    This version is not availble from GitHub or conda anymore; Keeping this version as the 1000G profiles from Brett Trost are based on this.

Scripts: 
R scripts: 
    1. are fetched from `crg` repo
    2. ~/crg/crg.ehdn.sh: runs the EHDN with many checks/searches for the relevant bam files based on run type: exome/genome and cre/crg/crg2 pipeline directory structures
    3. ~/crg/ehdn_report.sh: main script that uses all the R packages from ~/crg/str, 1000G EHDN profiles, and generates script
        a. ~/crg/str/DBSCAN.EHdn.parallel.R: library(dbscan), library(ggplot2), library(data.table), library(parallel), library(doSNOW)
        b. ~/crg/str//mergeExpansions.R: library(data.table), library(GenomicRanges), library(ggplot2), library(cowplot), library(Biostrings)
        c. 


Python scripts:
    ~/crg/str/compare_anchored_irrs.py: depends on "core" package local present
    ~/crg/str/find_outliers.py: numpy, matplotlib
    ~/crg/str/generate_EH_genotype_table.generic.py: argparse, glob, BTlib(local), collections, pandas, path
    ~/crg/str/combine_counts.py: argparse, json, core(local)
    ~/crg/str/BTlib.py: re, statistics, string, random, itertools, copy, collections, numpy, pandas, copy, docx
    ~/crg/str/add_gene+threshold_to_EH_column_headings2.py: argparse
    ~/crg/str/format_for_annovar.py: collections, re
    ~/crg/str/format_from_annovar.py: pandas, xlsxwriter
    


