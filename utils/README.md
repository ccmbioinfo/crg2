# Helpful scripts related to exome/genome processing

1. ```bam_slice_by_gene.py```: Extract a BAM file slice for a specific gene region or coordinates with configurable flanking sequence (source: https://github.com/naumenko-sa/bioscripts).
2. ```bam.coverage.base_wise.stat.py```: Calculate base-wise coverage stats for a BAM file (source: https://github.com/naumenko-sa/bioscripts). 
3. ```bam.coverage.less20.sh```: Report how many bases are below 20X coverage per exon in a BAM file (source: https://github.com/naumenko-sa/bioscripts).
4. ```bam.coverage.per_exon.sh```: Reports the distribution of lowest coverage in exons from a BAM file: (source: https://github.com/naumenko-sa/bioscripts).
5. ```check_dup_rate```: Calculate and report PCR duplication rates for all samples in a family.
6. ```copy_reports_to_minio_bcbio.sh```: Transfer bcbio pipeline analysis reports to MinIO storage (requires ccmmarvin credentials).
7. ```copy_reports_to_minio_DeepExomeSeq.sh```: Transfer deep exome sequencing reports to MinIO storage (requires ccmmarvin credentials).
8. ```copy_reports_to_minio.sh```: Transfer crg2 pipeline reports to MinIO storage (requires ccmmarvin credentials).
9. ```download_clinvar.sh```: Download and update ClinVar VCF files for both GRCh37 and GRCh38 genome builds. Add `0 0 15 * * sh ~/crg2/utils/download_clinvar.sh` to your crontab (`crontab -e`) to update ClinVar on the 15th of every month. 
10. ```gen_sample_sheet.py```: Create a sample sheet from TCAG genome data for use with crg2/parser_genomes.py.
11. ```mv2results.py```: Organize and move completed exome/genome analysis folders to the results directory.
12. ```wgs_gene_coverage/get_gene_coverage.sh```: Generate detailed coverage metrics for exons and introns of specified genes using mosdepth on WGS CRAM files.
