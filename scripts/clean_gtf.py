"""
Takes and cleans ensembl and refseq gtfs and gffs respectively (GRCh37), saving them as csvs.
wget -qO- http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz  | gunzip -c > Homo_sapiens.GRCh37.87.gtf
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz | gunzip -c > GRCh37_latest_genomic.gff
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_assembly_report.txt
Source: https://github.com/ccmbioinfo/sampletracker2stager/blob/master/variants/clean_gtf.py, Delvin So
"""

import pandas as pd
import argparse
import pyranges


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ensembl_gtf",
        type=str,
        help="unzipped ensembl gtf file",
    )
    parser.add_argument(
        "--refseq_gff3",
        type=str,
        help="unzipped refseq gff file",
    )
    parser.add_argument(
        "--refseq_assembly",
        type=str,
        help="refseq assembly report",
    )
    return parser


if __name__ == "__main__":

    args = get_parser().parse_args()

    # clean ensembl gtf
    ensembl = pyranges.read_gtf(args.ensembl_gtf)

    ensembl_ss = ensembl[(ensembl.Feature == "gene")]
    ensembl_ss = ensembl_ss[~(ensembl_ss.Chromosome.str.startswith("GL"))]

    ensembl_df = ensembl_ss.df[
        ["Chromosome", "Start", "End", "Source", "gene_id"]
    ].copy()
    ensembl_df.columns = map(str.lower, ensembl_df.columns)
    ensembl_df = ensembl_df.rename(columns={"gene_id": "name"})

    ensembl_df.to_csv("Homo_sapiens.GRCh38.110.gtf_subset.csv", index=False)

    # clean refseq gff (for some reason gtf throws errors with pyranges)
    refseq = pyranges.read_gff3(args.refseq_gff3)

    refseq_ss = refseq[(refseq.Feature == "gene")]
    refseq_df = refseq_ss.df[["Chromosome", "Start", "End", "Source", "Name"]].copy()
    refseq_df.columns = map(str.lower, refseq_df.columns)
    refseq_df.rename(columns={"chromosome": "chromosome_refseq"}, inplace=True)

    # map refseq chromosomes (e.g. accession NC_000001.10 refers to chr1)
    assembly = pd.read_csv(
        args.refseq_assembly,
        sep="\t",
        comment="#",
        names=[
            "Sequence-Name",
            "Sequence-Role",
            "Assigned-Molecule",
            "Assigned-Molecule-Location/Type",
            "GenBank-Accn",
            "Relationship",
            "RefSeq-Accn",
            "Assembly-Unit",
            "Sequence-Length",
            "UCSC-style-name",
        ],
    )
    assembly = assembly[["RefSeq-Accn", "Sequence-Name"]]
    chrom_dict = {}
    for accn, chr in zip(assembly["RefSeq-Accn"], assembly["Sequence-Name"]):
        chrom_dict[accn] = chr

    refseq_df["chromosome"] = [
        chrom_dict[accn] for accn in refseq_df["chromosome_refseq"].values
    ]

    # remove non-canonical chromosomes
    chrom = [str(i) for i in range(1, 23)] + ["MT", "X", "Y"]
    refseq_df = refseq_df[refseq_df["chromosome"].isin(chrom)]

    refseq_df = refseq_df.drop(columns=["chromosome_refseq"])
    refseq_df = refseq_df[["chromosome", "start", "end", "source", "name"]]

    refseq_df.to_csv("GRCh38_latest_genomic.gff_subset.csv", index=False)

    print("Done")
