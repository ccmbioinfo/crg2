1. ```concat_sample_tsv.sh```: combine coverage reports for multiple genes into one report for each sample. used when many genes are queried.
2. ```get_gene_coverage.sh```: generate wgs cram coverage report across exons/introns for given gene using mosdepth. 
3. ```get_gene_exon_intron.R```: generate BED with exon/intron coordiantes for query gene. used in get_gene_coverage.sh.
4. ```get_intron_gtf.sh```: get intron coordinates using bedtools. output used i
n parse_gtf.sh.
5. ```parse_gtf.R```: generate wgs cram coverage report across exons/introns for given gene using mosdepth.
6. ```slurm_hpf_gene_coverage.sh```: generate coverage reports given a gene list and cram directory.
