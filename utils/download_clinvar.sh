#!/bin/bash
# to run this script automatically on a monthly basis, add the following line to your crontab using `crontab -e`:
# 0 0 15 * * sh ~/crg2/utils/download_clinvar.sh
DATE=`date +%Y-%m-%d`
DIR_HG19=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/variation/
DIR_HG38=/hpf/largeprojects/ccmbio/nhanafi/c4r/downloads/databases/

# hg19
mv ${DIR_HG19}/clinvar.vcf.gz ${DIR_HG19}/clinvar_old/clinvar.${DATE}.vcf.gz
rm ${DIR_HG19}/clinvar.vcf.gz.tbi ${DIR_HG19}/clinvar.vcf.gz.md5

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz > ${DIR_HG19}/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi > ${DIR_HG19}/clinvar.vcf.gz.tbi
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.md5 > ${DIR_HG19}/clinvar.vcf.gz.md5
md5sum -c ${DIR_HG19}/clinvar.vcf.gz.md5

# hg38
mv ${DIR_HG38}/clinvar.vcf.gz ${DIR_HG38}/clinvar_old/clinvar.${DATE}.vcf.gz
rm ${DIR_HG38}/clinvar.vcf.gz.tbi ${DIR_HG38}/clinvar.vcf.gz.md5

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz > ${DIR_HG38}/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi > ${DIR_HG38}/clinvar.vcf.gz.tbi
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.md5 > ${DIR_HG38}/clinvar.vcf.gz.md5
md5sum -c ${DIR_HG38}/clinvar.vcf.gz.md5

