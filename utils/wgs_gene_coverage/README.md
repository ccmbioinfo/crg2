# Gene coverage scripts

1. ```get_gene_coverage.sh```: Generate cram coverage report from WES or WGS across exons/introns for given gene using mosdepth. 
2. ```get_gene_exon_intron.py```: Generate BED with exon/intron coordiantes for query gene. Used in get_gene_coverage.sh.
3. ```parse_gtf.R```: Generate BED with exon/intron coordiantes for the hg19 genome from RefSeq GTF file. Output is used in get_gene_exon_intron.R to get coordinates for specific genes.
4. ```intron_region_coverage_for_WES.py```: Calculate coverage at 10,20,50,100,1000 bp into intron sequences (useful for WES)
