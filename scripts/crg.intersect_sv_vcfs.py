import pandas as pd
from os import mkdir, path, environ
from shutil import move
from pathlib import Path
from datetime import date
import sys

from SVRecords import SVGrouper, SVAnnotator


def make_exon_gene_set(protein_coding_genes):
    df = pd.read_csv(protein_coding_genes, sep="\t")
    return set(df[df.columns[5]])


def filter_report(svreport):
    svreport = pd.read_csv(svreport, sep="\t", encoding="utf-8")
    filtered_report = svreport[svreport["NUM_SVTOOLS"] >= 2]
    return filtered_report


def main(
    reference,
    protein_coding_genes,
    exon_bed,
    hgmd_db,
    hpo,
    exac,
    omim,
    biomart,
    gnomad,
    sv_counts,
    outfile_name,
    vcfs,
    report_type,
):
    SVScore_cols = [
        "variants/SVLEN",
        "variants/SVSCORESUM",
        "variants/SVSCOREMAX",
        "variants/SVSCORETOP5",
        "variants/SVSCORETOP10",
        "variants/SVSCOREMEAN",
    ]
    MetaSV_col = "variants/NUM_SVTOOLS"
    HPO_cols = ["N_UNIQUE_HPO_TERMS", "HPO Features", "N_GENES_IN_HPO", "Genes in HPO"]
    Protein_coding_genes_col = "Protein-coding Ensembl Gene ID"
    protein_coding_ENSG = make_exon_gene_set(protein_coding_genes)

    print("Grouping like structural variants ...")
    try:
        sv_records = SVGrouper(
            vcfs, report_type=report_type, ann_fields=SVScore_cols + [MetaSV_col]
        )
    except KeyError:
        sv_records = SVGrouper(vcfs, report_type=report_type, ann_fields=SVScore_cols)
    sample_cols = [col for col in sv_records.df.columns if col != MetaSV_col]
    sample_genotype_cols = [col for col in sample_cols if col.endswith("_GENOTYPE")]

    print("Identifying protein coding genes ...")
    sv_records.df["Protein-coding Ensembl Gene ID"] = sv_records.df[
        "Ensembl Gene ID"
    ].apply(
        lambda gene_list: [gene for gene in gene_list if gene in protein_coding_ENSG]
    )

    print("Annotating structural variants ...")
    ann_records = SVAnnotator(reference, exon_bed, hgmd_db, hpo, exac, omim, biomart)
    sv_records.df = ann_records.annotate_genes(sv_records.df, Protein_coding_genes_col)

    for sv_count in sv_counts:
        prefix = Path(sv_count).stem
        sv_records.df = ann_records.annotate_counts(sv_count, sv_records, prefix=prefix)

    sv_records.df = ann_records.annotsv(reference,sv_records.df)
    sv_records.df = ann_records.calc_exons_spanned(sv_records.df, exon_bed)
    sv_records.df = ann_records.annotate_gnomad(gnomad, sv_records)
    sv_records.df = ann_records.annotate_hgmd(hgmd_db, sv_records.df)
    ann_records.add_decipher_link(sv_records.df)
    sv_records.df = ann_records.calc_exon_boundaries(
        sv_records.df.reset_index(), exon_bed
    )

    if not set(HPO_cols).issubset(set(sv_records.df.columns)):
        for col in HPO_cols:
            sv_records.df[col] = "na"

    # format and rearrange the columns
    if MetaSV_col not in sv_records.df.columns:
        sv_records.df["variants/NUM_SVTOOLS"] = "1"

    sv_report_cols = (
        [
            col
            for col in sample_cols
            if col
            not in set(
                SVScore_cols + sample_genotype_cols + ["N_SAMPLES", "Ensembl Gene ID"]
            )
        ]
        + [
            "variants/SVLEN",
        ]
        + [MetaSV_col, Protein_coding_genes_col]
        + [
            "BioMart Associated Gene Name",
            "EXONS_SPANNED",
            "Genes in HPO",
            "HPO Features",
            "Genes in OMIM",
            "OMIM Phenotypes",
            "OMIM Inheritance",
            "N_GENES_IN_HPO",
            "N_UNIQUE_HPO_TERMS",
            "N_GENES_IN_OMIM",
            "Canadian_MSSNG_parent_SVs.Manta.counts",
            "Canadian_MSSNG_parent_SVs.Manta.counts_SV",
            "Canadian_MSSNG_parent_SVs.LUMPY.counts",
            "Canadian_MSSNG_parent_SVs.LUMPY.counts_SV",
            "DGV_GAIN_IDs",
            "DGV_GAIN_n_samples_with_SV",
            "DGV_GAIN_n_samples_tested",
            "DGV_GAIN_Frequency",
            "DGV_LOSS_IDs",
            "DGV_LOSS_n_samples_with_SV",
            "DGV_LOSS_n_samples_tested",
            "DGV_LOSS_Frequency",
            "gnomAD_AF",
            "gnomAD_SV",
            "gnomAD_AN",
            "gnomAD_AC",
            "gnomAD_N_HOMREF",
            "gnomAD_N_HET",
            "gnomAD_N_HOMALT",
            "gnomAD_FREQ_HOMREF",
            "gnomAD_FREQ_HET",
            "gnomAD_FREQ_HOMALT",
            "gnomAD_GRPMAX_AF",
            "gnomAD_FILTER",
            "DDD_disease",
            "DDD_mode",
            "DDD_pmids",
            "Genes in HGMD",
            "HGMD disease",
            "HGMD descr",
            "HGMD JOURNAL_DETAILS",
            "ExAC syn_z",
            "ExAC mis_z",
            "ExAC lof_z",
            "ExAC pLI",
        ]
        + [col for col in SVScore_cols if col != "variants/SVLEN"]
        + [
            "nearestLeftExonBoundary",
            "nearestLeftExonDistance",
            "nearestRightExonBoundary",
            "nearestRightExonDistance",
            "DECIPHER_LINK",
        ]
    )
    # DGV CNVs not relevant for BNDs, nor are SVscores
    remove_from_bnd_report = (
        ["variants/SVLEN", MetaSV_col]
        + [
            "DGV_GAIN_IDs",
            "DGV_GAIN_n_samples_with_SV",
            "DGV_GAIN_n_samples_tested",
            "DGV_GAIN_Frequency",
        ]
        + [
            "DGV_LOSS_IDs",
            "DGV_LOSS_n_samples_with_SV",
            "DGV_LOSS_n_samples_tested",
            "DGV_LOSS_Frequency",
        ]
        + [col for col in SVScore_cols if col != "variants/SVLEN"]
    )

    bnd_report_cols = [
        col for col in sv_report_cols if col not in remove_from_bnd_report
    ]

    if report_type == "SV":
        sv_records.df = sv_records.df[sv_report_cols]
    else:
        sv_records.df = sv_records.df[bnd_report_cols]

    sv_records.df.columns = sv_records.df.columns.str.replace("variants/", "")
    sv_records.df = sv_records.df.drop_duplicates()

    # set missing values in numeric columns to 0, and '.' in non-numeric columns
    numeric = [
        "N_GENES_IN_HPO",
        "N_UNIQUE_HPO_TERMS",
        "N_GENES_IN_OMIM",
        "gnomAD_AF",
        "gnomAD_AN",
        "gnomAD_AC",
        "gnomAD_N_HOMREF",
        "gnomAD_N_HET",
        "gnomAD_N_HOMALT",
        "gnomAD_FREQ_HOMREF",
        "gnomAD_FREQ_HET",
        "gnomAD_FREQ_HOMALT",
        "gnomAD_GRPMAX_AF",
    ]
    non_numeric = [col for col in sv_records.df.columns if col not in numeric]
    for col in numeric:
        sv_records.df[col] = [
            val if val == val else "0" for val in sv_records.df[col].tolist()
        ]
    for col in non_numeric:
        print(col)
        sv_records.df[col] = [
            "." if val == "na" else val for val in sv_records.df[col].tolist()
        ]
        sv_records.df[col] = [
            "." if val == "" else val for val in sv_records.df[col].tolist()
        ]
        sv_records.df[col] = [
            "." if val == "nan" else val for val in sv_records.df[col].tolist()
        ]

    print("Writing results to file ...")
    sv_records.write(outfile_name)


if __name__ == "__main__":

    report_dir = snakemake.output[0]
    environ['TMPDIR'] = report_dir
    metasv = [vcf for vcf in snakemake.input if "metasv" in vcf]
    manta = [vcf for vcf in snakemake.input if "manta" in vcf]

    # create unfiltered SV report
    out_report = "{}.unfiltered.wgs.sv.v{}.{}.tsv".format(
        snakemake.params.project,
        snakemake.params.PIPELINE_VERSION,
        date.today().strftime("%Y-%m-%d"),
    )

    main(
        snakemake.params.reference,
        snakemake.params.protein_coding_genes,
        snakemake.params.exon_bed,
        snakemake.params.hgmd,
        snakemake.params.hpo,
        snakemake.params.exac,
        snakemake.params.omim,
        snakemake.params.biomart,
        snakemake.params.gnomad,
        [snakemake.params.mssng_manta_counts, snakemake.params.mssng_lumpy_counts],
        out_report,
        metasv,
        "SV",
    )

    # create report dir
    mkdir(report_dir)

    # create filtered SV report
    filtered_report = filter_report(out_report)
    out_filtered_report = "{}/{}.wgs.sv.v{}.{}.tsv".format(
        report_dir,
        snakemake.params.project,
        snakemake.params.PIPELINE_VERSION,
        date.today().strftime("%Y-%m-%d"),
    )
    filtered_report.to_csv(
        out_filtered_report, sep="\t", encoding="utf-8", na_rep=".", index=False
    )

    # create manta BND report
    bnd_report = "{}/{}.BND.wgs.sv.v{}.{}.tsv".format(
        report_dir,
        snakemake.params.project,
        snakemake.params.PIPELINE_VERSION,
        date.today().strftime("%Y-%m-%d"),
    )

    main(
        snakemake.params.reference,
        snakemake.params.protein_coding_genes,
        snakemake.params.exon_bed,
        snakemake.params.hgmd,
        snakemake.params.hpo,
        snakemake.params.exac,
        snakemake.params.omim,
        snakemake.params.biomart,
        snakemake.params.gnomad,
        [snakemake.params.mssng_manta_counts, snakemake.params.mssng_lumpy_counts],
        bnd_report,
        manta,
        "BND",
    )

    move(out_report, report_dir)
