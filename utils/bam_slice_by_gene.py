import argparse
import glob
import os
import requests
import pysam
import subprocess

# creates bam slices for a specified gene and assembly for all alignments in alignment path
# bam slices are output to /path/to/alignments/bam_slice/
# usage: python3 get_gene_coordinates.py -g <gene name> -f <flank size> -a <assembly, either hg19 or hg38> -p  </path/to/alignments>


def query_gene(gene, flank, assembly):
    print(f"Querying mygene.info for gene {gene} coordinates in assembly {assembly}")
    if assembly == "hg19":
        r = requests.get(
            f"https://mygene.info/v3/query?q=symbol:{gene}&species=human&field=genomic_pos_hg19"
        )
        coordinates = r.json()["hits"][0]["genomic_pos_hg19"]
    else:
        # hg38
        r = requests.get(
            f"https://mygene.info/v3/query?q=symbol:{gene}&species=human&field=genomic_pos"
        )
        coordinates = r.json()["hits"][0]["genomic_pos"]
    if type(coordinates) == list:
        coordinates = coordinates[0]
    chr, start, end = coordinates["chr"], coordinates["start"], coordinates["end"]
    flank = int(flank)
    start = start - flank
    end = end + flank
    coordinates = f"{chr}:{start}-{end}"
    print(f"{gene} coordinates: {coordinates}")

    return coordinates


def get_alignments(alignment_folder):
    alignments = glob.glob(f"{alignment_folder}/*ram")
    return alignments


def create_bam_slices(coordinates, alignment, ref, bam_slice_dir):
    sample = os.path.basename(alignment).split(".")[0]
    coordinates_name = coordinates.replace(":", "-")
    pysam.view(
        "-b",
        alignment,
        coordinates,
        "-T",
        ref,
        "-o",
        f"{bam_slice_dir}/{sample}_{coordinates_name}.bam",
        catch_stdout=False,
    )
    pysam.index(f"{bam_slice_dir}/{sample}_{coordinates_name}.bam")


if __name__ == "__main__":
    description = """Creates a bam slice for a specified gene
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-g",
        "--gene",
        type=str,
        required=True,
        help="Gene name (e.g. CDK2)",
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
        help="Integer that specifies the length of padding sequence around the gene",
    )

    parser.add_argument(
        "-p",
        "--path",
        type=str,
        required=True,
        help="Path to directory containing alignments (BAMs or CRAMs)",
    )
    args = parser.parse_args()
    gene = args.gene
    flank = args.flank
    assembly = args.assembly
    alignment_folder = args.path
    wd = os.getcwd()
    if "srv" in wd:
        ref = "/srv/shared/data/dccdipg/genomes/GRCh37d5/GRCh37d5.fa"
    else:
        ref = "/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa"

    bam_slice_dir = os.path.join(alignment_folder, "bam_slice")
    if not os.path.isdir(bam_slice_dir):
        os.mkdir(bam_slice_dir)

    # get gene coordinates
    coordinates = query_gene(gene, flank, assembly)
    alignments = get_alignments(alignment_folder)
    for a in alignments:
        create_bam_slices(coordinates, a, ref, bam_slice_dir)
