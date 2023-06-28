# Helpful scripts related to exome/genome processing

1. ```bam_slice_crg2.sh```: create a bam slice for a particular genomic region including specified flank.
2. ```check_dup_rate```: check duplication rate for each member of a family.
3. ```copy_reports_to_minio_bcbio.sh```: copy reports from bcbio runs to MinIO. Must be run as ccmmarvin.
4. ```copy_reports_to_minio_DeepExomeSeq.sh```: copy reports from deep exome sequencing runs to MinIO. Must be run as ccmmarvin.
5. ```copy_reports_to_minio.sh```: copy reports from crg2 runs to MinIO. Must be run as ccmmarvin.
6. ```gen_sample_sheet.py```: generate a sample sheet for a batch of genomes downloaded from TCAG for input to crg2/parser_genomes.py.
7. ```mv2results.py```: move exome or genome analysis folders to results directory.
8. ```cram_gene_coverage/get_gene_coverage.sh```: generate cram coverage report across exons/introns for given gene using mosdepth.
