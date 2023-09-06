import allel
import numpy as np
import pandas as pd
from pybedtools import BedTool
from SVRecords import SVAnnotator
from pathlib import Path
from datetime import date

# based on Dennis Kao's SVGrouper module: https://github.com/ccmbioinfo/crg/blob/master/SVRecords/SVGrouper.py

sv_counts = [snakemake.params.mssng_manta_counts]


def split_Ensembl_ids(id_list):
    new_list = []
    for id in id_list:
        if "-" in id:
            new_list.extend(id.split("-"))
        elif "&" in id:
            new_list.extend(id.split("&"))
        else:
            new_list.append(id)
    return new_list


def make_exon_gene_set(protein_coding_genes=snakemake.params.protein_coding_genes):
    df = pd.read_csv(protein_coding_genes, sep="\t")
    return set(df[df.columns[5]])


def make_ref_bedtool(sv_df):
    return BedTool(list(sv_df.index.values))


def annotate_counts(sv_df, counts, prefix="COUNT", reciprocal_overlap=0.5):
    print(
        "Annotating structural variants with those seen in %s based on a %f reciprocal overlap ..."
        % (counts, reciprocal_overlap)
    )

    cols = ["COUNT_CHROM", "COUNT_START", "COUNT_END", "COUNT_SVTYPE", "COUNT"]

    count_df = pd.read_csv(counts, sep="\t", dtype="str").astype(str)
    count_bed = BedTool(count_df.itertuples(index=False))

    sample_sv = make_ref_bedtool(sv_df)

    ann_df = sample_sv.intersect(
        count_bed, wa=True, wb=True, F=reciprocal_overlap, f=reciprocal_overlap
    ).to_dataframe(
        names=[
            "CHROM",
            "POS",
            "END",
            "SVTYPE",
        ]
        + cols
    )
    ann_df[
        [
            "CHROM",
            "POS",
            "END",
            "SVTYPE",
        ]
    ] = ann_df[
        [
            "CHROM",
            "POS",
            "END",
            "SVTYPE",
        ]
    ].astype(
        str
    )  # reference dataframe is typecasted as string
    ann_df = ann_df.drop(ann_df[ann_df["SVTYPE"] != ann_df["COUNT_SVTYPE"]].index)

    ann_df = ann_df.groupby(["CHROM", "POS", "END", "SVTYPE"]).agg(
        {
            "COUNT_CHROM": "first",
            "COUNT_SVTYPE": "first",
            "COUNT_START": "min",
            "COUNT_END": "max",
            "COUNT": "sum",
        }
    )
    ann_df["COUNT_SV"] = ann_df[["COUNT_CHROM", "COUNT_START", "COUNT_END"]].apply(
        lambda x: "{}:{}-{}".format(x[0], x[1], x[2]), axis=1
    )
    ann_df = ann_df.drop(
        columns=["COUNT_CHROM", "COUNT_START", "COUNT_END", "COUNT_SVTYPE"]
    )

    ann_df.columns = ann_df.columns.str.replace("COUNT", prefix)
    df = sv_df.join(ann_df)
    df[[prefix]] = df[[prefix]].fillna(0)

    return df


def annotate_gnomad(gnomad, sv_df):
    print("Annotating structural variants with those seen in gnomAD_SV")

    gnomad_cols = [
        "CHROM",
        "START",
        "END",
        "NAME",
        "SVTYPE",
        "AN",
        "AC",
        "AF",
        "N_HOMREF",
        "N_HET",
        "N_HOMALT",
        "FREQ_HOMREF",
        "FREQ_HET",
        "FREQ_HOMALT",
        "POPMAX_AF",
    ]
    gnomad_ann_cols = [
        "gnomAD_SVTYPE",
        "gnomAD_AN",
        "gnomAD_AC",
        "gnomAD_AF",
        "gnomAD_N_HOMREF",
        "gnomAD_N_HET",
        "gnomAD_N_HOMALT",
        "gnomAD_FREQ_HOMREF",
        "gnomAD_FREQ_HET",
        "gnomAD_FREQ_HOMALT",
        "gnomAD_POPMAX_AF",
    ]
    gnomad_df = pd.read_csv(snakemake.params.gnomad, sep="\t", dtype="str").astype(str)
    gnomad_df.columns = gnomad_df.columns.str.replace("#", "")
    gnomad_df.columns = gnomad_df.columns.str.strip()
    gnomad_df = gnomad_df[gnomad_cols]
    gnomad_bed = BedTool(gnomad_df.itertuples(index=False))

    sample_sv = make_ref_bedtool(sv_df)

    ann_df = (
        sample_sv.intersect(
            gnomad_bed,
            wa=True,
            wb=True,
        )
        .to_dataframe(
            names=[
                "CHROM",
                "POS",
                "END",
                "SVTYPE",
                "gnomAD_CHROM",
                "gnomAD_START",
                "gnomAD_END",
                "gnomAD_ID",
            ]
            + gnomad_ann_cols
        )
        .astype(str)
    )
    ann_df = ann_df.drop(ann_df[ann_df["SVTYPE"] != ann_df["gnomAD_SVTYPE"]].index)

    ann_df["gnomAD_SV"] = ann_df[["gnomAD_CHROM", "gnomAD_START", "gnomAD_END"]].apply(
        lambda x: "{}:{}-{}".format(x[0], x[1], x[2]), axis=1
    )
    ann_df = ann_df.drop(columns=["gnomAD_CHROM", "gnomAD_START", "gnomAD_END"])
    ann_df = ann_df.groupby(["CHROM", "POS", "END", "SVTYPE"]).agg(list)
    ann_df = ann_df[ann_df.columns].applymap(lambda cell: " & ".join(cell))

    return sv_df.join(ann_df)


def main():
    vcf_dict = allel.read_vcf(
        snakemake.input.vcf[0],
        ["*"],  # extract all fields
        numbers={"ANN": 1000},  # ?
        transformers=allel.ANNTransformer(),  # post-processing of data from SNPEF
    )

    protein_coding_ENSG = make_exon_gene_set(snakemake.params.protein_coding_genes)
    Protein_coding_genes_col = "Protein-coding Ensembl Gene ID"

    # grab alt allele
    vcf_dict["variants/ALT"] = [alt[0] for alt in vcf_dict["variants/ALT"]]
    vcf_dict["variants/ANN_Annotation"] = [
        list(filter(None, list(set(impact))))
        for impact in vcf_dict["variants/ANN_Annotation"]
    ]
    sample_sv_fields = (
        [
            "variants/CHROM",
            "variants/POS",
            "variants/END",
            "variants/SVTYPE",
            "variants/SVTYPE",
        ]
        + [
            "variants/ALT",
            "variants/ANN_Annotation",
            "calldata/AD",
            "calldata/DP",
        ]
        + ["calldata/GT", "variants/ANN_Gene_ID", "samples"]
    )

    parse_fields = list(set(sample_sv_fields))

    # drop un-needed fields from vcf, cannot pass in parse_fields to read_vcf() because ANN_gene_id is unknown until ANNTransformer runs
    for key in list(vcf_dict.keys()):
        if key not in parse_fields:
            vcf_dict.pop(key)

    # remove empty strings, split on delimited characters, then join using comma
    # by default, specifying numbers=1000 creates 1000 elements, with most being empty
    vcf_dict["variants/ANN_Gene_ID"] = [
        list(filter(None, ann)) for ann in vcf_dict["variants/ANN_Gene_ID"]
    ]
    vcf_dict["variants/ANN_Gene_ID"] = [
        split_Ensembl_ids(id_list)
        if any("&" in id for id in id_list) or any("-" in id for id in id_list)
        else id_list
        for id_list in vcf_dict["variants/ANN_Gene_ID"]
    ]
    vcf_dict["variants/ANN_Gene_ID"] = [
        ",".join(list(set(id_list))) if isinstance(id_list, list) else id_list
        for id_list in vcf_dict["variants/ANN_Gene_ID"]
    ]

    # parse genotype fields
    gts = {}
    names = vcf_dict["samples"]
    for index, name in enumerate(names):
        gts[str(index)] = []

    for gt in vcf_dict["calldata/GT"]:
        for index, g in enumerate(gt):
            if set(g) == {0}:
                g = "-"
            elif set(g) == {1}:
                g = "hom"
            else:
                g = "het"
            gts[str(index)].append(g)

    for index, sample in enumerate(vcf_dict["samples"]):
        vcf_dict[f"{sample}_GENOTYPE"] = gts[str(index)]

    # parse AD field
    ads = {}
    names = vcf_dict["samples"]
    for index, name in enumerate(names):
        ads[str(index)] = []

    for ad in vcf_dict["calldata/AD"]:
        for index, ad in enumerate(ad):
            ad = ad[0]
            ads[str(index)].append(ad)

    for index, sample in enumerate(vcf_dict["samples"]):
        vcf_dict[f"{sample}_AD"] = ads[str(index)]

    # parse DP field
    dps = {}
    names = vcf_dict["samples"]
    for index, name in enumerate(names):
        dps[str(index)] = []

    for dp in vcf_dict["calldata/DP"]:
        for index, dp in enumerate(dp):
            dps[str(index)].append(dp)

    for index, sample in enumerate(vcf_dict["samples"]):
        vcf_dict[f"{sample}_DP"] = dps[str(index)]

    # remove unneccessary fields
    for key in "samples", "calldata/GT", "calldata/AD", "calldata/DP":
        del vcf_dict[key]

    df = pd.DataFrame(vcf_dict)

    df = df.rename(
        {
            "variants/ALT": "ALT",
            "variants/ANN_Gene_ID": "Ensembl Gene ID",
            "variants/CHROM": "CHROM",
            "variants/POS": "POS",
            "variants/SVTYPE": "INSTYPE",
            "variants/ANN_Annotation": "Impact",
        },
        axis=1,
    )
    df["END"] = df["POS"]
    df["SVTYPE"] = "INS"
    df["POS"] = df["POS"].astype(str)
    df["END"] = df["END"].astype(str)

    df = df.set_index(["CHROM", "POS", "END", "SVTYPE"])

    print("Identifying protein coding genes ...")
    genes = []
    for g in df["Ensembl Gene ID"]:
        try:
            g_list = g.split(",")
        except:
            print(g)
            g_list.append("NA")
        genes.append(g_list)
    df["variants/ANN_Gene_ID"] = genes
    df["Protein-coding Ensembl Gene ID"] = df["variants/ANN_Gene_ID"].apply(
        lambda gene_list: [gene for gene in gene_list if gene in protein_coding_ENSG]
    )

    print("Annotating structural variants against genes, hgmd, hpo, exac, and omim")
    HPO_cols = ["N_UNIQUE_HPO_TERMS", "HPO Features", "N_GENES_IN_HPO", "Genes in HPO"]
    ann_records = SVAnnotator(
        snakemake.params.exon_bed,
        snakemake.params.hgmd_db,
        snakemake.params.hpo,
        snakemake.params.exac,
        snakemake.params.omim,
        snakemake.params.biomart,
    )
    df = ann_records.annotate_genes(df, Protein_coding_genes_col)

    # annotate against MSSNG Manta calls (not Lumpy as Lumpy set doesn't include INS)
    for sv_count in sv_counts:
        prefix = Path(sv_count).stem
        df = annotate_counts(df, sv_count, prefix=prefix)

    # calculate number of exons spanned
    df = ann_records.calc_exons_spanned(df, snakemake.params.exon_bed)

    # annotate against gnomAD SV
    df = annotate_gnomad(snakemake.params.gnomad, df)

    # annotate against hgmd
    df = ann_records.annotate_hgmd(snakemake.params.hgmd_db, df)

    # add DECIPHER URL
    ann_records.add_decipher_link(df)

    # calculate distance to exon boundaries
    df = ann_records.calc_exon_boundaries(df.reset_index(), snakemake.params.exon_bed)

    # if HPO terms not provided, make these columns "na"
    if not set(HPO_cols).issubset(set(df.columns)):
        for col in HPO_cols:
            df[col] = "na"

    # parse impacts
    df["Impact"] = df["Impact"].apply(lambda impact: ("|").join(impact))

    # subset columns
    sample_cols = [
        col for col in df.columns if "GENOTYPE" in col or "_AD" in col or "_DP" in col
    ]
    sv_report_cols = (
        [col for col in sample_cols]
        + [Protein_coding_genes_col]
        + [
            "INSTYPE",
            "Impact",
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
            "gnomAD_POPMAX_AF",
            "Genes in HGMD",
            "HGMD disease",
            "HGMD descr",
            "HGMD JOURNAL_DETAILS",
            "ExAC syn_z",
            "ExAC mis_z",
            "ExAC lof_z",
            "ExAC pLI",
        ]
        + [
            "nearestLeftExonBoundary",
            "nearestLeftExonDistance",
            "nearestRightExonBoundary",
            "nearestRightExonDistance",
            "DECIPHER_LINK",
        ]
    )

    df = df[sv_report_cols]

    # parse datatypes
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
        "gnomAD_POPMAX_AF",
    ]
    non_numeric = [col for col in df.columns if col not in numeric]
    for col in numeric:
        print(col)
        df[col] = [val if val == val else "0" for val in df[col].tolist()]
    for col in non_numeric:
        print(col)
        df[col] = ["." if val == "na" else val for val in df[col].tolist()]
        df[col] = ["." if val == "" else val for val in df[col].tolist()]
        df[col] = ["." if val == "nan" else val for val in df[col].tolist()]

    # get report file name
    report_name = "{}.wgs.MELT.{}.csv".format(
        snakemake.params.family, date.today().strftime("%Y-%m-%d")
    )
    # write report
    df.to_csv("report/MELT/" + report_name)


main()
