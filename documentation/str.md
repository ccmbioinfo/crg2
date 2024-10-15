# ExpansionHunterDenovo version: 
    EHDN v0.7.0 is not available from GitHub or conda anymore, and the 1000G profiles from Brett Trost are based on this, so keeping this version; tool path in 'config.yaml' is updated to
    `config['tools']['ehdn']: "/hpf/largeprojects/ccm_dccforge/dccdipg/Common/crg2-non-conda-tools/EHDN.TCAG/ExpansionHunterDenovo-v0.7.0"`

# Description of scripts used before the changes:
## BASH and R scripts:

1. scripts are fetched from `crg/str` repo
2. ~/crg/crg.ehdn.sh: runs  EHDN with many checks/searches for the relevant bam files based on run type: exome/genome and cre/crg/crg2 pipeline directory structures
3. ~/crg/ehdn_report.sh: main script that uses all the R packages from ~/crg/str, 1000G EHDN profiles, and generates script
    - ~/crg/str/DBSCAN.EHdn.parallel.R: library(dbscan), library(ggplot2), library(data.table), library(parallel), library(doSNOW)
    - ~/crg/str//mergeExpansions.R: library(data.table), library(GenomicRanges), library(ggplot2), library(cowplot), library(Biostrings)
    

## Python scripts:
1. ~/crg/str/compare_anchored_irrs.py: depends on "core" package  present locally
2. ~/crg/str/find_outliers.py: numpy, matplotlib
3. ~/crg/str/generate_EH_genotype_table.generic.py: argparse, glob, BTlib(local), collections, pandas, path
4. ~/crg/str/combine_counts.py: argparse, json, core(local)
5. ~/crg/str/BTlib.py: re, statistics, string, random, itertools, copy, collections, numpy, pandas, copy, docx
6. ~/crg/str/add_gene+threshold_to_EH_column_headings2.py: argparse
7. ~/crg/str/format_for_annovar.py: collections, re
8. ~/crg/str/format_from_annovar.py: pandas, xlsxwriter


# Summary of changes related to issue 142:
1. Moved all required scripts from `~/crg` and `~/cre` to `~/crg2/scripts/str`
2. Removed hard-coded paths from str related scripts and added them to `config.yaml` to be passed as arguments via snakemake params.
3. Moved following annotations and tools to /hpf/largeprojects/ccm_dccforge/dccdipg/Common/:
    annotation/ExpansionHunterDenovo/UCSC_simple_repeats_hg19_coord_motif.tsv
    crg2-non-conda-tools/EHDN.TCAG/ExpansionHunterDenovo-v0.7.0
    annotation/ExpansionHunterDenovo/1000G_JSON
    updated annotation/ExpansionHunterDenovo/manifest.1000G.txt
4. Split the single "EHdn_report" rule into four, with separate env for each rule 
    - EHDN_mark_outliers: 
    - EHDN_DBSCAN_outlier
    - EHDN_merge_expansions
    - EHDN_annovar
4. Created `envs/ehdn-dbscan.yaml` for rules "EHDN_DBSCAN_outlier" and "EHDN_merge_expansions". The version of the R-packages are not the same as in `ccmmarvin`. If I constrain them to versions available in 'ccmmarvin`, then conda fails to create env due to conflicts. The versions installed on ccmmarvin were done manually in 2020, not via conda, 
5. Fixed minor bugs, and removed unwanted yaml. 
6. Edited above R scripts to remove hard-coded file paths, add command-line arguments, and removed  suffixing output with date, as the Snakemake rule requires the output names be known before execution (using dynamic only works if all outputs from a rule are dynamic)
7. Tested with NA12878



