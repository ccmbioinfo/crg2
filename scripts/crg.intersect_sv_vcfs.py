import argparse
import pandas as pd
from pathlib import Path
from SVRecords import SVGrouper, SVAnnotator

def make_exon_gene_set(protein_coding_genes):
    df = pd.read_csv(protein_coding_genes, sep="\t")
    return(set(df[df.columns[5]]))

def main(protein_coding_genes, exon_bed, hgmd_db, hpo, exac, omim, biomart, gnomad, sv_counts, outfile_name, vcfs):
    SVScore_cols = ['variants/SVLEN', 'variants/SVSCORESUM', 'variants/SVSCOREMAX', 'variants/SVSCORETOP5', 'variants/SVSCORETOP10', 'variants/SVSCOREMEAN',]
    MetaSV_col = 'variants/NUM_SVTOOLS'
    HPO_cols = [ "N_UNIQUE_HPO_TERMS", "HPO Features", "N_GENES_IN_HPO", "Genes in HPO" ]
    Protein_coding_genes_col = "Protein-coding Ensembl Gene ID"
    protein_coding_ENSG = make_exon_gene_set(protein_coding_genes)

    print("Grouping like structural variants ...")
    sv_records = SVGrouper(vcfs, ann_fields=SVScore_cols + [MetaSV_col])
    sample_cols = [ col for col in sv_records.df.columns if col != MetaSV_col ]
    sample_genotype_cols = [col for col in sample_cols if col.endswith('_GENOTYPE')]

    print("Identifying protein coding genes ...")
    sv_records.df['Protein-coding Ensembl Gene ID'] = sv_records.df['Ensembl Gene ID'].apply(lambda gene_list: [gene for gene in gene_list if gene in protein_coding_ENSG])

    print('Annotating structural variants ...')
    ann_records = SVAnnotator(exon_bed, hgmd_db, hpo, exac, omim, biomart)
    sv_records.df = ann_records.annotate_genes(sv_records.df, Protein_coding_genes_col)

    for sv_count in sv_counts:
        prefix = Path(sv_count).stem
        sv_records.df = ann_records.annotate_counts(sv_count, sv_records, prefix=prefix)
    
    sv_records.df = ann_records.annotsv(sv_records.df)
    sv_records.df = ann_records.calc_exons_spanned(sv_records.df, exon_bed)
    sv_records.df = ann_records.annotate_gnomad(gnomad, sv_records)
    sv_records.df = ann_records.annotate_hgmd(hgmd_db, sv_records.df)
    ann_records.add_decipher_link(sv_records.df)
    sv_records.df = ann_records.calc_exon_boundaries(sv_records.df.reset_index(), exon_bed)

    if not set(HPO_cols).issubset(set(sv_records.df.columns)):
        for col in HPO_cols:
            sv_records.df[col] = "na"

    # format and rearrange the columns
    sv_records.df = sv_records.df[ [col for col in sample_cols if col not in set(SVScore_cols + sample_genotype_cols + ['N_SAMPLES', 'Ensembl Gene ID']) ] + \
    [ 'variants/SVLEN', ] + \
    [ MetaSV_col ] + \
    [ Protein_coding_genes_col ] + \
    [ "BioMart Associated Gene Name", "EXONS_SPANNED", ] + \
    [ "Genes in HPO", "HPO Features", ] + \
    [ "Genes in OMIM", "OMIM Phenotypes", "OMIM Inheritance", ] + \
    [ "N_GENES_IN_HPO", "N_UNIQUE_HPO_TERMS", "N_GENES_IN_OMIM", ] + \
    [ "Canadian_MSSNG_parent_SVs.Manta.counts", "Canadian_MSSNG_parent_SVs.Manta.counts_SV", "Canadian_MSSNG_parent_SVs.LUMPY.counts", "Canadian_MSSNG_parent_SVs.LUMPY.counts_SV"] + \
    [ "DGV_GAIN_IDs", "DGV_GAIN_n_samples_with_SV", "DGV_GAIN_n_samples_tested", "DGV_GAIN_Frequency", ] + \
    [ "DGV_LOSS_IDs", "DGV_LOSS_n_samples_with_SV", "DGV_LOSS_n_samples_tested", "DGV_LOSS_Frequency", ] + \
    [ "gnomAD_AF", "gnomAD_SV", "gnomAD_AN", "gnomAD_AC", "gnomAD_N_HOMREF", "gnomAD_N_HET", "gnomAD_N_HOMALT", "gnomAD_FREQ_HOMREF", "gnomAD_FREQ_HET", "gnomAD_FREQ_HOMALT", "gnomAD_POPMAX_AF" ] + \
    [ "DDD_disease", "DDD_mode", "DDD_pmids", ] + \
    [ "Genes in HGMD", "HGMD disease", "HGMD descr", "HGMD JOURNAL_DETAILS" ] + \
    [ "ExAC syn_z", "ExAC mis_z", "ExAC lof_z", "ExAC pLI" ] + \
    [ col for col in SVScore_cols if col != 'variants/SVLEN'] + \
    [ "nearestLeftExonBoundary", "nearestLeftExonDistance", "nearestRightExonBoundary", "nearestRightExonDistance", ] + \
    [ "DECIPHER_LINK" ] ]
    sv_records.df.columns = sv_records.df.columns.str.replace('variants/','')
    sv_records.df = sv_records.df.drop_duplicates()

    # [ "DDD_SV", "DDD_DUP_n_samples_with_SV", "DDD_DUP_Frequency", "DDD_DEL_n_samples_with_SV", "DDD_DEL_Frequency" ]

    # set missing values in numeric columns to 0, and '.' in non-numeric columns
    numeric =  [ "N_GENES_IN_HPO", "N_UNIQUE_HPO_TERMS", "N_GENES_IN_OMIM","gnomAD_AF", "gnomAD_AN", "gnomAD_AC", "gnomAD_N_HOMREF", "gnomAD_N_HET", "gnomAD_N_HOMALT", "gnomAD_FREQ_HOMREF", "gnomAD_FREQ_HET", "gnomAD_FREQ_HOMALT", "gnomAD_POPMAX_AF" ]
    non_numeric = [col for col in sv_records.df.columns if col not in numeric]
    for col in numeric:
        sv_records.df[col] = [val if val == val else '0' for val in sv_records.df[col].tolist()]
    for col in non_numeric:
        sv_records.df[col] = ['.' if val == 'na' else val for val in sv_records.df[col].tolist()]
        sv_records.df[col] = ['.' if val == '' else val for val in sv_records.df[col].tolist()]
        sv_records.df[col] = ['.' if val == 'nan' else val for val in sv_records.df[col].tolist()]


    print('Writing results to file ...')
    sv_records.write(outfile_name)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generates a structural variant report in CSV format for clincal research')
    parser.add_argument('-i', type=str, nargs='+', help='VCF files containing structural variants, e.g. -i 180_230.vcf 180_231.vcf', required=True)
    parser.add_argument('-protein_coding_genes', help='BED file containing protein coding gene regions with their ENSG id\'s', required=True)
    parser.add_argument('-exon_bed', help='BED file containing fixed exon positions', required=True)
    parser.add_argument('-hgmd', help='HGMD Pro database file', required=True, type=str)
    parser.add_argument('-hpo', help='Tab delimited file containing gene names and HPO terms', type=str)
    parser.add_argument('-exac', help='ExAC tab delimited file containing gene names and scores', type=str, required=True)
    parser.add_argument('-omim', help='OMIM tab delimited file containing gene names and scores', type=str, required=True)
    parser.add_argument('-biomart', help='TSV file from BiomaRt containing Ensemble gene ID, transcript ID, gene name, MIM gene id, HGNC id, EntrezGene ID', type=str, required=True)
    parser.add_argument('-gnomad', help='BED file from Gnomad containing structural variant coordinates and frequencies across populations', type=str, required=True)
    parser.add_argument('-sv_counts', nargs='+', help='List of BED files containing structural variants and their frequencies. Can be used to annotate with various populations and variant callers', required=False)
    parser.add_argument('-overlap', help='Recipricol overlap to group a structural variant by', type=float, default=0.5)
    parser.add_argument('-o', help='Output file name e.g. -o 180.sv.family.tsv', required=True, type=str)
    args = parser.parse_args()

    if len(args.i) == 0:
        ValueError('Please enter the path to some vcf\'s following the -i flag')
    else:
        main(args.protein_coding_genes, args.exon_bed, args.hgmd, args.hpo, args.exac, args.omim, args.biomart, args.gnomad, args.sv_counts, args.o, args.i)
