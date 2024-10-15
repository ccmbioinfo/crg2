import pandas as pd, sys
import xlsxwriter

'''
format the outputs from ANNOVAR table_annovar.pl with OMIM and GNOMAD
merged.rare.EHdn.annovar.omim.hg19_multianno.txt
merged.rare.EHdn.annovar.gnomad.hg19_multianno.txt
'''

gnomad_file = sys.argv[1]
omim_file = sys.argv[2]
outfile = sys.argv[3]


#read OMIM annotated file as table and fix header
omim = pd.read_csv(omim_file, sep="\t", header=[0,1], index_col=None)
omim.columns = [j if i.startswith("Otherinfo") else i for i, j in omim.columns ]
omim.rename(columns=lambda x: x.replace(".refGene",""), inplace=True)
omim_col = ['omim_inheritance', 'omim_phenotype' ]

#read gnoMAD annotated file as table and fix header
gnomad = pd.read_csv(gnomad_file, sep="\t", header=[0,1], index_col=None)
gnomad.columns = [j if i.startswith("Otherinfo") else i for i, j in gnomad.columns ]
gnomad.rename(columns=lambda x: x.replace(".refGene",""), inplace=True)

#sample-wise outlier and size column name for reporting
sample_outliers, size_outliers = [], []
for i in list(gnomad.columns):
    if i.startswith("outlier.") or i.startswith("1000G_outlier_count"):
        sample_outliers.append(i)
    if i.startswith("size."):
        size_outliers.append(i)
print(sample_outliers, size_outliers)
print(gnomad.columns)

gnomad_col = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func', 'Gene', 'motif' ] + sample_outliers + size_outliers + ['key', '1000G_freq_perc', 'GeneDetail','ExonicFunc', 'AAChange', 'transcript', 'oe_lof', 'oe_lof_upper', 'oe_mis', 'oe_mis_upper', 'pLI', 'pRec']
print(gnomad_col)
annot_rare = pd.merge(gnomad[gnomad_col],omim[omim_col],right_index=True, left_index=True)
print(annot_rare.head(5))

output_columns = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Function', 'Gene', 'motif', ]
output_columns += sample_outliers + size_outliers
output_columns += ['chr#start#end', 'a1000g_freq_perc', 'GeneDetail','ExonicFunc', 'AAChange', 'transcript', 'oe_lof', 'oe_lof_upper', 'oe_mis', 'oe_mis_upper', 'pLI', 'pRec', 'omim_inheritance', 'omim_phenotype' ]
if not ".xlsx" in outfile:
    outfile = outfile + ".xlsx"
annot_rare.to_excel(outfile, index=False)
