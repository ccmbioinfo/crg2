import pandas as pd
import numpy as np
from pysam import VariantFile
import argparse
from collections import defaultdict
from pybedtools import BedTool
from sigfig import round
from datetime import date


def rename_SV_cols(annotsv_df):
    annotsv_df.rename(
        columns={
            "SV_chrom": "CHROM",
            "SV_start": "POS",
            "SV_end": "END",
            "SV_length": "SVLEN",
            "SV_type": "SVTYPE",
        },
        inplace=True,
    )
    return annotsv_df


def annotsv_df_to_bed(annotsv_df):
    annotsv_df = annotsv_df[["CHROM", "POS", "END", "SVTYPE", "ID"]]
    bed = BedTool.from_dataframe(annotsv_df)
    return bed


def filter_length(row):
    length = row.SVLEN
    if pd.isna(length):
        return True
    elif abs(length) >= 50:
        return True
    else:
        return False


def apply_filter_length(annotsv_df):
    annotsv_df_filtered = annotsv_df.apply(lambda row: filter_length(row), axis=1)
    return annotsv_df_filtered


def filter_benign(row):
    svtype = row.SVTYPE
    loss_AF_max = row.B_loss_AFmax
    gain_AF_max = row.B_gain_AFmax
    ins_AF_max = row.B_ins_AFmax
    inv_AF_max = row.B_inv_AFmax
    if svtype == "DEL":
        if loss_AF_max == "nan":
            return True
        else:
            return False
    elif svtype == "DUP":
        if gain_AF_max == "nan":
            return True
        else:
            return False
    elif svtype == "INS":
        if ins_AF_max == "nan":
            return True
        else:
            return False
    elif svtype == "INV":
        if inv_AF_max == "nan":
            return True
        else:
            return False
    # no benign annotations for BNDs
    else:
        return True


def apply_filter_benign(annotsv_df):
    annotsv_df_filtered = annotsv_df.apply(lambda row: filter_benign(row), axis=1)
    return annotsv_df_filtered


def merge_full_split_annos(annotsv_df):
    annotsv_split = annotsv_df[annotsv_df["Annotation_mode"] == "split"]
    annotsv_full = annotsv_df[annotsv_df["Annotation_mode"] == "full"]
    DDD_cols = ["DDD_status", "DDD_mode", "DDD_consequence", "DDD_disease", "DDD_pmid"]
    GenCC_cols = ["GenCC_disease", "GenCC_moi", "GenCC_classification", "GenCC_pmid"]
    transcript_cols = [
        "Tx",
        "Tx_start",
        "Tx_end",
        "Overlapped_tx_length",
        "Overlapped_CDS_length",
        "Overlapped_CDS_percent",
        "Frameshift",
        "Exon_count",
        "Location",
        "Location2",
        "Dist_nearest_SS",
        "Nearest_SS_type",
        "Intersect_start",
        "Intersect_end",
    ]
    cols_of_interest = DDD_cols + GenCC_cols + transcript_cols
    annotsv_split_agg = (
        annotsv_split.groupby("AnnotSV_ID")
        .agg(
            {
                "GenCC_disease": ";".join,
                "GenCC_moi": ";".join,
                "GenCC_classification": ";".join,
                "GenCC_pmid": ";".join,
                "Tx": ";".join,
                "Tx_start": ";".join,
                "Tx_end": ";".join,
                "Overlapped_tx_length": ";".join,
                "Overlapped_CDS_length": ";".join,
                "Overlapped_CDS_percent": ";".join,
                "Frameshift": ";".join,
                "Exon_count": ";".join,
                "Location": ";".join,
                "Location2": ";".join,
                "Nearest_SS_type": ";".join,
                "Dist_nearest_SS": ";".join,
                "Intersect_start": ";".join,
                "Intersect_end": ";".join,
                "DDD_status": ";".join,
                "DDD_mode": ";".join,
                "DDD_consequence": ";".join,
                "DDD_disease": ";".join,
                "DDD_pmid": ";".join,
            }
        )
        .reset_index()
    )

    annotsv_full = annotsv_full[
        [col for col in annotsv_full.columns if col not in cols_of_interest]
    ]
    annotsv_merge = annotsv_split_agg.merge(
        annotsv_full, how="outer", on="AnnotSV_ID"
    ).reset_index()
    return annotsv_merge


def add_hpo(hpo, gene):
    try:
        gene = gene.split(";")
    except AttributeError:
        return "NA"
    terms = []
    for g in gene:
        try:
            term = str(hpo[hpo["Gene symbol"] == g]["Features"].values[0])
            term = term.replace("; ", ";").split(";")
            term = list(set(term))
            for t in term:
                terms.append(t)
        except IndexError:
            pass
    if len(terms) == 0:
        return "nan"
    else:
        terms = ",".join(terms)
        return terms


def add_omim(omim_df, gene):
    gene = gene.split(";")
    phenos = []
    inheritance = []
    for g in gene:
        try:
            phenos.append(
                str(omim_df[omim_df["gene_name"] == g]["omim_phenotype"].values[0])
            )
            inheritance.append(
                str(omim_df[omim_df["gene_name"] == g]["omim_inheritance"].values[0])
            )
        except IndexError:
            pass
    phenos = ",".join(phenos)
    inheritance = ",".join(inheritance)
    return [phenos, inheritance]


def get_genotype(sample_GT_AD_DP):
    genotype = sample_GT_AD_DP.split(":")[0]
    if genotype == "0/1":
        return "het"
    elif genotype == "1/1":
        return "hom"
    elif genotype == "./.":
        return genotype
    elif genotype == "0/0":
        return "-"
    else:
        return genotype


def get_alt_depth(sample_GT_AD_DP):
    alt_depth = sample_GT_AD_DP.split(":")[1].split(",")[1]
    return alt_depth


def get_depth(sample_GT_AD_DP):
    depth = sample_GT_AD_DP.split(":")[2]
    return depth


def get_exon_counts(annotsv_df, exon_bed):
    # original functions from SVRecords/SVAnnotator.py
    # this function includes the PB SV ID, because there can be multiple SVs with the same CHROM POS END SVTYPE that have different ALTs
    # these are uniquely identified by the PB SV ID
    exon_counts = defaultdict(
        int
    )  # incrementing dict and setting value in df is faster than incrementing values in the df
    exon_ref = BedTool(exon_bed)
    sample_bedtool = BedTool(
        list(annotsv_df.reset_index()[["CHROM", "POS", "END", "SVTYPE", "ID"]].values)
    )

    for interval in sample_bedtool.intersect(exon_ref, wa=True):
        exon_counts[
            (
                str(interval.chrom),
                str(interval.start),
                str(interval.stop),
                str(interval[3]),
                str(interval[4]),
            )
        ] += 1

    count_df = pd.Series(exon_counts).to_frame()
    count_df = count_df.reset_index()
    count_df.columns = ["CHROM", "POS", "END", "SVTYPE", "ID", "EXONS_SPANNED"]
    count_df = count_df.astype(str)
    annotsv_df = pd.merge(
        annotsv_df, count_df, how="left", on=["CHROM", "POS", "END", "SVTYPE", "ID"]
    ).fillna(value={"EXONS_SPANNED": 0})
    return annotsv_df


def annotate_pop_svs(annotsv_df, pop_svs, cols):
    # intersect annotsv and population SV bed file
    annotsv_bed = annotsv_df_to_bed(annotsv_df)
    pop_svs = pd.read_csv(pop_svs, sep="\t")
    pop_svs_bed = BedTool.from_dataframe(pop_svs)
    intersect = annotsv_bed.intersect(
        pop_svs_bed, wa=True, wb=True, F=0.5, f=0.5
    ).to_dataframe()
    intersect.columns = [
        "CHROM",
        "POS",
        "END",
        "SVTYPE",
        "ID",
        "CHROM_pop",
        "POS_pop",
        "END_pop",
        "SVTYPE_pop",
    ] + cols
    # popSV and sample SV must be same type
    intersect = intersect[intersect["SVTYPE"] == intersect["SVTYPE_pop"]]
    pop_name = cols[0].split("_")[0]
    # e.g. make a column with SV details, e.g DEL:1:25266309-25324509
    intersect[f"{pop_name}_SV"] = intersect[
        ["CHROM_pop", "POS_pop", "END_pop", "SVTYPE_pop"]
    ].apply(lambda x: f"{x[3]}:{x[0]}:{x[1]}-{x[2]}", axis=1)
    cols.append(f"{pop_name}_SV")
    intersect = intersect[["CHROM", "POS", "END", "SVTYPE", "ID"] + cols]
    # round AFs
    try:
        AF_col = [col for col in cols if "AF" in col][0]
        intersect[AF_col] = [
            round(float(af), sigfigs=3) for af in intersect[AF_col].values
        ]
    except IndexError:
        pass
    # group by SV, joining annotation columns
    intersect = intersect.astype(str)
    intersect = (
        intersect.groupby(["CHROM", "POS", "END", "SVTYPE", "ID"])[cols]
        .agg({col: "; ".join for col in cols})
        .reset_index()
    )
    # get max allele frequency
    try:
        AF_col = [col for col in cols if "AF" in col][0]
        intersect[f"{pop_name}_maxAF"] = intersect[AF_col].apply(
            lambda x: max([float(af) for af in x.split("; ")])
        )
        cols.append(f"{pop_name}_maxAF")
    # get max allele counts for C4R
    except IndexError:
        count_cols = ["C4R_AC", "seen_in_C4R_count"]
        for col in count_cols:
            intersect[f"{col}_max"] = intersect[col].apply(
                lambda x: max([int(ac) for ac in x.split("; ")])
            )
        cols.append(f"{col}_max")

    # merge population AF dataframe with annotSV df
    annotsv_pop_svs = pd.merge(
        annotsv_df,
        intersect,
        how="left",
        on=["CHROM", "POS", "END", "SVTYPE", "ID"],
    ).fillna(value={col: 0 for col in cols})
    return annotsv_pop_svs


def main(df, omim, hpo, vcf, prefix, exon_bed, cmh, hprc, gnomad, inhouse):
    # filter out SVs < 50bp
    df_len = df[apply_filter_length(df)]
    df_len = df_len.astype(str)
    # merge full and split AnnotSV annos
    df_merge = merge_full_split_annos(df_len)
    sample_cols = [col for col in df.columns if "PB" in col]
    # extract genotype and alt allele depth
    for sample in sample_cols:
        df_merge[f"{sample}_GT"] = [
            get_genotype(row[sample]) for index, row in df_merge.iterrows()
        ]
        df_merge[f"{sample}_AD"] = [
            get_alt_depth(row[sample]) for index, row in df_merge.iterrows()
        ]
        df_merge[f"{sample}_DP"] = [
            get_depth(row[sample]) for index, row in df_merge.iterrows()
        ]
    gt_cols = [col for col in df_merge.columns if "_GT" in col]
    ad_cols = [col for col in df_merge.columns if "_AD" in col]
    dp_cols = [col for col in df_merge.columns if "_DP" in col]

    # filter out benign SVs with AF > 1%
    # df_merge_notbenign = df_merge[apply_filter_benign(df_merge)]

    # add HPO terms by gene matching
    df_merge["HPO"] = [add_hpo(hpo, gene) for gene in df_merge["Gene_name"].values]

    # add OMIM phenos and inheritance by gene matching
    df_merge["omim_phenotype"] = [
        add_omim(omim, gene)[0] for gene in df_merge["Gene_name"].values
    ]
    df_merge["omim_inheritance"] = [
        add_omim(omim, gene)[1] for gene in df_merge["Gene_name"].values
    ]

    # add SVs called by pbsv from Children's Mercy Hospital PacBio HiFi data
    cmh_cols = ["cmh_AF", "cmh_AC"]
    df_merge = annotate_pop_svs(df_merge, cmh, cmh_cols)

    # add SVs called by pbsv from HPRC PacBio HiFi data
    hprc_cols = ["hprc_AF", "hprc_AC", "hprc_HOM"]
    df_merge = annotate_pop_svs(df_merge, hprc, hprc_cols)

    # add gnomAD SVs
    gnomad_cols = [
        "gnomad_NAME",
        "gnomad_POPMAX_AF",
        "gnomad_AC",
        "gnomad_HOM",
    ]
    df_merge = annotate_pop_svs(df_merge, gnomad, gnomad_cols)

    # add C4R inhouse db SV counts
    inhouse_cols = [
        "C4R_ID",
        "C4R_REF",
        "C4R_ALT",
        "C4R_AC",
        "seen_in_C4R",
        "seen_in_C4R_count",
    ]
    df_merge = annotate_pop_svs(df_merge, inhouse, inhouse_cols)

    # add exon counts
    df_merge = get_exon_counts(df_merge, exon_bed)
    print(df_merge.columns)

    # define columns to be included in report
    report_columns = (
        [
            "CHROM",
            "POS",
            "END",
            "SVLEN",
            "SVTYPE",
            "ID",
            "INFO",
            "Gene_name",
            "Gene_count",
            "omim_phenotype",
            "omim_inheritance",
            "DDD_mode",
            "DDD_consequence",
            "DDD_disease",
            "HPO",
        ]
        + gt_cols
        + ad_cols
        + dp_cols
        + ["Tx", "Frameshift", "EXONS_SPANNED", "Nearest_SS_type", "Dist_nearest_SS"]
        + cmh_cols
        + hprc_cols
        + gnomad_cols
        + inhouse_cols
        + [
            "DDD_HI_percent",
            "ExAC_delZ",
            "ExAC_dupZ",
            "ExAC_cnvZ",
            "ExAC_synZ",
            "ExAC_misZ",
            "ExAC_pLI",
            "CytoBand",
            "RE_gene",
            "TAD_coordinate",
            "ENCODE_experiment",
            "Repeat_type_left",
            "Repeat_type_right",
            "SegDup_left",
            "SegDup_right",
            "ENCODE_blacklist_characteristics_left",
            "ENCODE_blacklist_characteristics_right",
        ]
    )

    df_merge = df_merge[report_columns]
    df_merge = df_merge.drop(columns=["C4R_REF", "C4R_ALT"])
    df_merge["Gene_name"] = df_merge["Gene_name"].replace("nan", ".")
    df_merge["omim_phenotype"].fillna("nan", inplace=True)
    df_merge["omim_inheritance"].fillna("nan", inplace=True)
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    df_merge.to_csv(f"{prefix}.rare.hpo.{today}.csv", index=False)

    # make a dictionary of variants with SVLEN >=50 and BNDs from vcf
    # validate that the only SVs missing from annotSV merged and SV >=50 bp are on non-canonical chr
    vcf_variant_dict = {}
    i = 0
    for record in vcf.fetch():
        i += 1
        try:
            length = abs(float(record.info["SVLEN"][0]))
            if length >= 50:
                vcf_variant_dict[record.id] = length
        except KeyError:
            if record.info["SVTYPE"] == "BND":
                vcf_variant_dict[record.id] = 0

    vcf_ID = vcf_variant_dict.keys()
    annotSV_ID = list(set(df_merge["ID"]))
    print("SVs in pbsv that are not present in AnnotSV tsv:")
    print(vcf_ID - annotSV_ID)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generates a structural variant report using an AnnotSV-annotated pbsv VCF generated from PacBio HiFi WGS"
    )
    parser.add_argument("-annotsv", type=str, help="AnnotSV tsv file", required=True)
    parser.add_argument(
        "-hpo",
        help="Tab delimited file containing gene names and HPO terms",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-omim",
        help="OMIM tab delimited file containing gene names and scores",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-vcf",
        help="pbsv vcf file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-exon",
        help="exon bed file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-cmh",
        help="Children's Mercy Hospital SVs in bed format",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-hprc",
        help="HPRC SVs in bed format",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-gnomad",
        help="gnomad SVs in bed format",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-inhouse",
        help="C4R inhouse database",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    df = pd.read_csv(args.annotsv, sep="\t", low_memory=False)
    df = rename_SV_cols(df)
    prefix = args.annotsv.strip(".tsv")
    vcf = VariantFile(args.vcf)
    omim = pd.read_csv(args.omim, sep="\t")
    omim = omim[pd.notnull(omim["omim_phenotype"])]
    hpo = pd.read_csv(
        args.hpo,
        comment="#",
        skip_blank_lines=True,
        sep="\t",
        encoding="ISO-8859-1",
        engine="python",
    )
    # Phenotips TSV has a space in column name: " Gene symbol"
    hpo.columns = hpo.columns.str.strip()
    main(
        df,
        omim,
        hpo,
        vcf,
        prefix,
        args.exon,
        args.cmh,
        args.hprc,
        args.gnomad,
        args.inhouse,
    )
