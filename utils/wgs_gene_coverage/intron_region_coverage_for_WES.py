import pandas as pd
import statistics as stats
import sys
import tabix 

# usage: intron_coverage.py <gene_name> <mosdepth_per-base_cov> <gene_regions.bed>
# where mosdepth_per-base_cov is the per-base coverage bed file (eg sample.per-base.bed.gz) output by mosdepth when '-n' argument is excluded
# gene_regions.bed is a bed file of intron and exon coordinates for transcripts of a specified gene, generated by crg2/utils/wgs_gene_coverage/get_gene_coverage.sh 
# the script outputs the mean coverage for the start and end regions of each intron (e.g. 10bp into the intron and 10bp at the end of the intron)

def get_mean_cov(cov_tbx, chr, start, end, region):
    # get average coverage at first x bases (specificed by region argument) and last x bases of an intron
    chr = str(chr)
    region = int(region)
    r1_end = start + region
    r2_start = end - region
    cov = []
    for q in cov_tbx.query(chr, start, r1_end):
        cov.append(int(q[3]))
    for q in cov_tbx.query(chr, r2_start, end):
        cov.append(int(q[3]))
    median_cov = stats.mean(cov)
    
    return median_cov

# read in command-line coverage
gene = sys.argv[1]
cov_by_base = sys.argv[2]
gene_coord = sys.argv[3]
gene_coord = pd.read_csv(gene_coord, sep="\t")

# exclude introns
gene_coord_intron = gene_coord[gene_coord["TYPE"] == "intron"]

# use tabix to query coveraged BED file
cov_tbx = tabix.open(cov_by_base)

for region in ["10", "20", "50", "100", "500"]:
    col_name = f"region_{region}"
    print(col_name)
    gene_coord_intron[col_name] = gene_coord_intron.apply(lambda row: get_mean_cov(cov_tbx, row['CHR'], row[f'START'], row[f'END'], region), axis=1)


gene_coord_intron.to_csv(f"gene_coverage/{gene}_intron_region_coverage.tsv", index=False, sep="\t")

