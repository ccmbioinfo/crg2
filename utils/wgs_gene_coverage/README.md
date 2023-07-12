# WGS gene coverage scripts

1. ```get_gene_coverage.sh```: Generate wgs cram coverage report across exons/introns for given gene using mosdepth. 
2. ```get_gene_exon_intron.R```: Generate BED with exon/intron coordiantes for query gene. Used in get_gene_coverage.sh.
3. ```parse_gtf.R```: Generate BED with exon/intron coordiantes for the hg19 genome from RefSeq GTF file. Output is used in get_gene_exon_intron.R to get coordinates for specific genes.
