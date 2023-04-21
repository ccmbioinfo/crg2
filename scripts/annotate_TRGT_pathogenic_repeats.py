import allel
import argparse
import numpy as np
import pandas as pd
from pysam import VariantFile
from collections import defaultdict
import os


def parse_sample_field(vcf_dict, vcf_df, fieldname, fieldtype):
    field_dict = {}

    for record in vcf_dict.fetch():
        for sample in record.samples:
            value = record.samples[sample][fieldname]
            value = [str(v) for v in value]
            name = record.samples[sample].name
            name = name.split(".")[0]
            if name not in field_dict.keys():
                field_dict[name] = [value]
            else:
                field_dict[name].append(value)
    for sample in record.samples:
        name = record.samples[sample].name
        name = name.split(".")[0]
        vcf_df[f"{name}_{fieldtype}"] = [
            ("|").join(field) for field in field_dict[name]
        ]

    return vcf_df


def recode_genes(disease_thresholds):
    disease_thresholds.loc[disease_thresholds["Gene"] == "DMPK/DMPKas", "Gene"] = "DMPK"
    disease_thresholds.loc[
        disease_thresholds["Gene"] == "ATXN8/ATXN8OS", "Gene"
    ] = "ATXN8"
    # disease_thresholds.loc[disease_thresholds["Gene"] == "NOTCH2NLC", "Gene"] = "NOTCH2NL"
    # in the TRGT vcf, one locus is annotated as NOTCH2NL, one as NOTCH2NLC- are they both NOTCH2NLC?
    disease_thresholds.loc[disease_thresholds["Gene"] == "CBL2", "Gene"] = "CBL"
    disease_thresholds.loc[disease_thresholds["Gene"] == "XYLT", "Gene"] = "XYLT1"
    disease_thresholds.loc[disease_thresholds["Gene"] == "TK2/BEAN", "Gene"] = "BEAN1"
    disease_thresholds.loc[disease_thresholds["Gene"] == "MARCH6", "Gene"] = "MARCHF6"
    disease_thresholds.loc[disease_thresholds["Gene"] == "C9orf72", "Gene"] = "C9ORF72"
    disease_thresholds.loc[disease_thresholds["Gene"] == "FMR1/FMR4", "Gene"] = "FMR1"
    disease_thresholds.loc[
        disease_thresholds["Gene"] == "LOC642361/NUTM2B-AS1", "Gene"
    ] = "NUTM2B-AS1"
    disease_thresholds.loc[disease_thresholds["Gene"] == "ZNF9/CNBP", "Gene"] = "CNBP"
    return disease_thresholds


def is_disease(motif_count, threshold, motif, allele_length):
    if "_" in motif_count:
        # ATNX8
        # motif with submotifs
        # calculate motif count: divide allele length by motif length
        if threshold == ">110":
            motif_len = len(motif)
            if allele_length == "None":
                motif_count = "."
            else:
                alleles = allele_length.split("|")
                a1, a2 = float(int(alleles[0]) / motif_len), float(
                    int(alleles[1]) / motif_len
                )
                motif_count = [a1, a2]
        else:
            # for other genes with submotifs, don't have reliable way to interpret pathogencity
            return None
    else:
        try:
            motif_count = [int(m) for m in motif_count.split("|")]
        except:
            motif_count = "."

    is_disease = False
    if threshold != threshold or motif_count == ".":
        # threshold is NaN
        is_disease = None
    elif "(" in threshold:
        is_disease = None
    elif "-" in threshold:
        try:
            min_thresh = int(threshold.split("-")[0])
        except ValueError:
            min_thresh = int(threshold.replace(">", "").split("-")[0])
        max_thresh = int(threshold.split("-")[1])
        for count in motif_count:
            if count >= min_thresh and count <= max_thresh:
                is_disease = True
    elif ">" in threshold:
        threshold = int(threshold.replace(">", ""))
        for count in motif_count:
            if count > int(threshold):
                is_disease = True
        for count in motif_count:
            if count > int(threshold):
                is_disease = True
    elif "," in threshold:
        path_counts = [int(t) for t in threshold.split(",")]
        for count in motif_count:
            if count in path_counts:
                is_disease = True
    else:
        for count in motif_count:
            if count >= int(threshold):
                is_disease = True
    return is_disease


def main(vcf, disease_thresholds):
    vcf_dict = allel.read_vcf(
        vcf,
        ["*"],
    )

    disease_thresholds = pd.read_csv(
        disease_thresholds,
        sep="\t",
    )

    vcf_dict["variants/ALT"] = ["|".join(alt) for alt in vcf_dict["variants/ALT"]]

    remove_keys = [
        "calldata/AP",
        "calldata/AL",
        "calldata/GT",
        "variants/FILTER_PASS",
        "variants/ID",
        "variants/altlen",
        "variants/is_snp",
        "calldata/ALCI",
        "calldata/AM",
        "calldata/MC",
        "calldata/MS",
        "calldata/SD",
        "samples",
    ]

    for key in remove_keys:
        vcf_dict.pop(key)

    vcf_df = pd.DataFrame(vcf_dict)

    # allel doesn't extract alt and ref support properly, so use pysam VariantFile instead
    vcf_dict_pysam = VariantFile(
        vcf,
    )

    vcf_df = parse_sample_field(vcf_dict_pysam, vcf_df, "AL", "allele_length")
    vcf_df = parse_sample_field(vcf_dict_pysam, vcf_df, "ALCI", "allele_CI")
    vcf_df = parse_sample_field(vcf_dict_pysam, vcf_df, "MC", "motif_count")
    vcf_df = parse_sample_field(vcf_dict_pysam, vcf_df, "MS", "motif_span")
    vcf_df = parse_sample_field(vcf_dict_pysam, vcf_df, "AM", "avg_methylation")
    vcf_df = parse_sample_field(vcf_dict_pysam, vcf_df, "SD", "spanning_read_support")
    vcf_df.columns = [col.replace("variants/", "") for col in vcf_df.columns]

    # recode genes names in disease threshold file for compatibility with TRGT gene names
    disease_thresholds = recode_genes(disease_thresholds)

    # extract disease threshold and disorder for each gene
    thresh_dict = {}
    disorder_dict = {}
    for index, row in disease_thresholds.iterrows():
        gene = row["Gene"]
        threshold = row["Disease threshold"]
        disorder = row["Disorder"]
        thresh_dict[gene] = threshold
        disorder_dict[gene] = disorder

    vcf_df["DISEASE_THRESHOLD"] = vcf_df["TRID"].map(thresh_dict)
    vcf_df["DISORDER"] = vcf_df["TRID"].map(disorder_dict)

    # fill in missing disease thresholds (from PRichmond's thesis)
    vcf_df.loc[vcf_df["TRID"] == "C11ORF80", "DISEASE_THRESHOLD"] = "500"
    vcf_df.loc[vcf_df["TRID"] == "AFF2", "DISEASE_THRESHOLD"] = "200"
    vcf_df.loc[vcf_df["TRID"] == "TMEM185A", "DISEASE_THRESHOLD"] = "300"

    # fill in missing disorders (OMIM)
    vcf_df.loc[
        vcf_df["TRID"] == "C11ORF80", "DISORDER"
    ] = "Spastic paraplegia 6, autosomal dominant"
    vcf_df.loc[
        vcf_df["TRID"] == "NIPA1", "DISORDER"
    ] = "Spastic paraplegia 6, autosomal dominant	"
    vcf_df.loc[
        vcf_df["TRID"] == "AFF2", "DISORDER"
    ] = "Intellectual developmental disorder, X-linked 109"

    # prep report for export
    allele_length_cols = [col for col in vcf_df.columns if "allele_length" in col]
    allele_CI_cols = [col for col in vcf_df.columns if "allele_CI" in col]
    spanning_cols = [col for col in vcf_df.columns if "spanning" in col]
    motif_count_cols = [col for col in vcf_df.columns if "motif_count" in col]
    methylation_cols = [col for col in vcf_df.columns if "methylation" in col]

    # add disease outlier column
    for col in motif_count_cols:
        vcf_df[f"DISEASE_PREDICTION_{col}"] = vcf_df.apply(
            lambda row: is_disease(
                row[col],
                row.DISEASE_THRESHOLD,
                row.MOTIFS,
                row[col.replace("motif_count", "allele_length")],
            ),
            axis=1,
        )
    disease_pred_cols = [col for col in vcf_df.columns if "PREDICTION" in col]
    vcf_df["DISEASE_PREDICTION"] = vcf_df[disease_pred_cols].apply(
        lambda row: "|".join(row.values.astype(str)), axis=1
    )
    for col in disease_pred_cols:
        vcf_df.drop(col, axis=1)

    report_cols = (
        [
            "CHROM",
            "POS",
            "END",
            "REF",
            "ALT",
            "TRID",
            "DISORDER",
            "MOTIFS",
            "STRUC",
            "DISEASE_THRESHOLD",
            "DISEASE_PREDICTION",
        ]
        + motif_count_cols
        + allele_length_cols
        + allele_CI_cols
        + spanning_cols
        + methylation_cols
    )

    vcf_df = vcf_df[report_cols]
    vcf_df = vcf_df.rename(columns={"TRID": "GENE"})
    outfile = vcf.replace(".vcf", "_TRGT_pathogenic_repeats.csv")
    vcf_df.to_csv(f"{outfile}", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generates a structural variant report using an AnnotSV-annotated pbsv VCF generated from PacBio HiFi WGS"
    )
    parser.add_argument("-vcf", type=str, help="TRGT VCF", required=True)

    args = parser.parse_args()
    vcf = args.vcf
    disease_thresholds = "/hpf/largeprojects/ccm_dccforge/dccdipg/Common/annotation/ExpansionHunter/tandem_repeat_disease_loci_v1.1.tsv"

    main(
        vcf,
        disease_thresholds,
    )
