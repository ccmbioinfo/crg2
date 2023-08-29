import pandas as pd
import os
import sys


def read_genome_BED():
    '''
    Read genome intron/exon coordinates (from parse_gft.R)
    '''
    genome_bed = pd.read_csv(
        "/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/hg19_introns_exons.bed",
        sep=" ",
        skiprows=1,
        low_memory=False,
        names=["seqid", "start", "end", "type", "gene_id", "transcript_id"],
    )
    return genome_bed


def get_valid_genes(df):
    '''
    Get list of unique genes from df
    '''
    unique_genes = df["gene_id"].unique()
    return unique_genes


def get_gene_BED(gene_name, bed):
    '''
    Get bed dataframe filtered for gene_name
    '''
    
    print(f"Input gene: {gene_name}")
    
    ## check if gene_name is valid
    valid_genes = get_valid_genes(bed)

    if gene_name in valid_genes:
        ## filter genome_bed for gene
        gene_bed = bed.loc[bed["gene_id"] == gene_name, :]
        print(f"Generating BED for: {gene_name}")
        return gene_bed

    ## print error if invalid gene_name
    else:
        print(f"'{gene_name}' is an invalid gene name!")
        sys.exit()


def concat_annotations(bed):
    '''
    Concatenate annotation columns in bed into "id" column (mosdepth only allows one metadata column)
    '''
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


def write_BED(bed, gene_name, project_dir):
    '''
    Make gene_coverage folder in project_dir and write bed file to folder
    '''
    target_dir = project_dir + "/gene_coverage/"

    if not os.path.exists(target_dir):
        # Create the directory
        os.mkdir(target_dir)
        print(f"Directory '{target_dir}' created.")
    else:
        print(f"Directory '{target_dir}' already exists.")

    bed.to_csv(f"{target_dir}/{gene_name}_pos.bed", sep="\t", header=False, index=False)



def main():
    ## Get variables from the command line
    gene_name=sys.argv[1]
    project_dir=sys.argv[2]
    
    ## Call functions
    genome_bed = read_genome_BED()
    gene_bed = get_gene_BED(gene_name=gene_name, bed=genome_bed)
    gene_bed_concat = concat_annotations(bed=gene_bed)
    write_BED(bed=gene_bed_concat, gene_name=gene_name, project_dir=project_dir)
    
    print(f"{gene_name} BED file saved to {project_dir}/gene_coverage/")


main()