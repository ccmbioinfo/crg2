import pandas as pd
import os
import sys
import jupyter_black

jupyter_black.load()

## Read genome intron/exon coordinates (from parse_gft.R)
def Read_Genome_BED():
    genome_bed = pd.read_csv(
        "/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/hg19_introns_exons.bed",
        sep=" ",
        skiprows=1,
        low_memory=False,
        names=["seqid", "start", "end", "type", "gene_id", "transcript_id"],
    )
    return genome_bed


## Get valid genes list
def Get_Valid_Genes(df):
    unique_genes = df["gene_id"].unique()
    return unique_genes


## Get BED dataframe filtered for a query gene
def Get_Gene_BED(gene_name, bed):
    print("Input gene: " + gene_name)

    ## check if gene_name is valid
    valid_genes = Get_Valid_Genes(bed)

    if gene_name in valid_genes:
        ## filter genome_bed for gene
        gene_bed = bed.loc[bed["gene_id"] == gene_name, :]
        print("Generating BED for: " + gene_name)
        return gene_bed

    ## print error if invalid gene_name
    else:
        print(f"'{gene_name}' is an invalid gene name!")
        sys.exit()


## Concatenate annotation columns into "id" column (mosdepth only allows one metadata column)
def Concat_Annotations(bed):
    bed_concat = bed.copy()

    bed_concat["id"] = (
        bed_concat["transcript_id"]
        + ";"
        + bed_concat["gene_id"]
        + ";"
        + bed_concat["type"]
    )

    bed_concat = bed_concat.drop(columns=["transcript_id", "gene_id", "type"])

    return bed_concat


## Make gene_coverage folder in project_dir and write BED file to folder
def Write_BED(bed, gene_name, project_dir):
    target_dir = project_dir + "/gene_coverage/"

    if not os.path.exists(target_dir):
        # Create the directory
        os.mkdir(target_dir)
        print(f"Directory '{target_dir}' created.")
    else:
        print(f"Directory '{target_dir}' already exists.")

    bed.to_csv(f"{target_dir}/{gene_name}_pos.bed", sep="\t", header=False, index=False)



def main():
    genome_bed = Read_Genome_BED()
    gene_bed = Get_Gene_BED(gene_name=sys.argv[1], bed=genome_bed)
    gene_bed_concat = Concat_Annotations(bed=gene_bed)

    Write_BED(bed=gene_bed_concat, gene_name=sys.argv[1], project_dir=sys.argv[2])
    print(sys.argv[1] + "BED file saved to " + sys.argv[2] + "/gene_coverage/")


main()