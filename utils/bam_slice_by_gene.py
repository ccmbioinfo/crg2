import argparse
import glob
import os
import requests
import pysam
import subprocess
import sys

# creates bam slices for a specified gene and assembly for all alignments in alignment path
# bam slices are output to <current workding directory/bam_slice/
# assumes you are running script from folder containing bams from which to create slices
# usage: python3 bam_slice_by_gene.py -g <gene name> -f <flank size> -a <assembly, either hg19 or hg38> -p  <platform>


def query_gene(gene, flank, assembly):
    # query mygene.info with gene name to retrieve gene coordinates
    print(f"Querying mygene.info for gene {gene} coordinates in assembly {assembly}")
    if assembly == "hg19":
        r = requests.get(
            f"https://mygene.info/v3/query?q=symbol:{gene}&species=human&field=genomic_pos_hg19"
        )
        print(r.json()["hits"])
        try:
            coordinates = r.json()["hits"][0]["genomic_pos_hg19"]
        except KeyError:
            coordinates = r.json()["hits"][1]["genomic_pos_hg19"]
        except IndexError:
            print(f"{gene} not found, exiting")
            sys.exit()
    else:
        # hg38
        r = requests.get(
            f"https://mygene.info/v3/query?q=symbol:{gene}&species=human&field=genomic_pos"
        )
        coordinates = r.json()["hits"][0]["genomic_pos"]
    if type(coordinates) != dict:
        for i in coordinates:
            chromosomes = [str(chr) for chr in range(1, 23)] + ["X", "Y"]
            if i["chr"] in chromosomes:
                coordinates = i
    chr, start, end = coordinates["chr"], coordinates["start"], coordinates["end"]
    flank = int(flank)
    start = start - flank
    end = end + flank
    coordinates = f"{chr}:{start}-{end}"
    if assembly == "hg38":
        coordinates = f"chr{coordinates}"
    print(f"{gene} coordinates: {coordinates}")

    return coordinates


def get_alignments(platform):
    # get bams or crams present in working directory
    if platform == "PacBio":
        patterns = ["*haplotagged.bam", "*asm*bam"] # we only want bam slices from haplotagged bam and assembled bam
        alignments = []
        for pattern in patterns:
            alignments.extend(glob.glob(pattern))
    else:
        alignments = glob.glob(f"*am")
    return alignments


def create_bam_slices(coordinates, alignment, ref, platform, gene=None):
    # use samtools to generate bam slices
    sample = os.path.basename(alignment).split(".")[0]
    if platform == "PacBio":
        if "asm" in alignment:
            sample = sample + "_asm"
        elif "haplotagged" in alignment:
            sample = sample + "_haplotagged"
    coordinates_name = coordinates.replace(":", "-")
    if gene:
        coordinates_name = f"{coordinates_name}-{gene}"
    pysam.view(
        "-b",
        alignment,
        coordinates,
        "-T",
        ref,
        "-o",
        f"bam_slice/{sample}_{coordinates_name}.bam",
        catch_stdout=False,
    )
    pysam.index(f"bam_slice/{sample}_{coordinates_name}.bam")


if __name__ == "__main__":
    description = """Creates a bam slice for a specified gene or coordinates provided
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-g",
        "--gene",
        type=str,
        required=False,
        help="Gene name (e.g. CDK2)",
    )
    parser.add_argument(
        "-c",
        "--coordinates",
        type=str,
        required=False,
        help="Coordinates, e.g. 14:35186001-35985000. Only used if gene name not provided. ",
    )
    parser.add_argument(
        "-f",
        "--flank",
        type=int,
        required=True,
        help="Integer that specifies the length of padding sequence around the gene",
    )
    parser.add_argument(
        "-a",
        "--assembly",
        type=str,
        choices=["hg19", "hg38"],
        required=True,
        help="Reference assembly version",
    )
    parser.add_argument(
        "-p",
        "--platform",
        type=str,
        required=True,
        choices=["Illumina", "PacBio"],
        help="Sequencing platform (Illumina or PacBio)"
    )

    args = parser.parse_args()
    gene = args.gene
    flank = args.flank
    assembly = args.assembly
    platform = args.platform

    if not gene:
        coordinates = args.coordinates.replace(",", "")
        if assembly == "hg19":
            coordinates = coordinates.replace("chr", "")
        elif assembly == "hg38":
            chr = coordinates.split(":")[0]
            if "chr" not in chr:
                coordinates = "chr" + coordinates
        chr = coordinates.split(":")[0]
        start = int(coordinates.split(":")[1].split("-")[0]) - flank
        end = int(coordinates.split(":")[1].split("-")[1]) + flank
        coordinates = f"{chr}:{start}-{end}"

    wd = os.getcwd()
    if "srv" in wd:
        ref = "/srv/shared/data/dccdipg/genomes/GRCh37d5/GRCh37d5.fa"
    else:
        ref = "/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa"

    if not os.path.isdir("bam_slice"):
        os.mkdir("bam_slice")

    alignments = get_alignments(platform)

    if gene:
        coordinates = query_gene(gene, flank, assembly)
        for a in alignments:
            create_bam_slices(coordinates, a, ref, platform, gene=gene)

    else:
        for a in alignments:
            create_bam_slices(coordinates, a, platform, ref)
