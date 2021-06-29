#!/usr/bin/env python3

import pandas as pd
import sys
from os.path import basename, splitext, join

#Usage: python3 add_hpo_terms_to_wes.py <phenotips_hpo_terms> <wes.report.csv> <outdir>

HPO_DF = pd.read_csv(sys.argv[1], comment='#', skip_blank_lines=True,\
	sep="\t",  encoding="ISO-8859-1", engine='python').drop(columns=["HPO IDs"])

#Phenotips TSV has a space in column name: " Gene symbol"
HPO_DF.columns = HPO_DF.columns.str.strip()
HPO_DF = HPO_DF.rename(columns={'Gene symbol': 'Gene Symbol'})

WES_REPORT = pd.read_csv(sys.argv[2], encoding="ISO-8859-1").set_index("Position")

OUT = join(sys.argv[3],"%s.w_hpoterms.tsv" % splitext(basename(sys.argv[2]))[0])

if "Ensembl_gene_id" in WES_REPORT.columns:
        HPO_DF = HPO_DF.set_index("Gene ID").drop(columns=["Gene Symbol"])
        OUT_DF = WES_REPORT.join(HPO_DF, on="Ensembl_gene_id")\
                .rename(columns={"Number of occurrences": "HPO_count", "Features": "HPO_terms"}).fillna("NA")
else:
        HPO_DF = HPO_DF.set_index("Gene Symbol").drop(columns=["Gene ID"])
       	OUT_DF = WES_REPORT.join(HPO_DF, on="Gene")\
                .rename(columns={"Number of occurrences": "HPO_count", "Features": "HPO_terms"}).fillna("NA")

OUT_DF.to_csv(OUT, sep="\t", encoding='utf-8')
