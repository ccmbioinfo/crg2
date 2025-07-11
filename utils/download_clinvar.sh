#!/bin/bash
# to run this script automatically on a monthly basis, add the following line to your crontab using `crontab -e`:
# 0 0 15 * * sh ~/crg2/utils/download_clinvar.sh
set -e
server=$1
DATE=`date +%Y-%m-%d`

if [ -z "$server" ]; then
    echo "Usage: $0 <server>"
    echo "server: hpf or g4rd"
    exit 1
fi

if [ "$server" == "hpf" ]; then
    DIR_HG19=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/variation/
    DIR_HG38=/hpf/largeprojects/ccmbio/nhanafi/c4r/downloads/databases/
elif [ "$server" == "g4rd" ]; then
    DIR_HG19=/srv/shared/data/dccdipg/annotation/vcfanno/
fi

# hg19
mv ${DIR_HG19}/clinvar.vcf.gz ${DIR_HG19}/clinvar_old/clinvar.${DATE}.vcf.gz
rm ${DIR_HG19}/clinvar.vcf.gz.tbi ${DIR_HG19}/clinvar.vcf.gz.md5

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz -O ${DIR_HG19}/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi -O ${DIR_HG19}/clinvar.vcf.gz.tbi
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.md5 -O ${DIR_HG19}/clinvar.vcf.gz.md5
md5=`awk '{print $1}'  ${DIR_HG19}/clinvar.vcf.gz.md5`
echo -e "${md5}\t${DIR_HG19}/clinvar.vcf.gz" > ${DIR_HG19}/clinvar.vcf.gz.md5
md5sum -c ${DIR_HG19}/clinvar.vcf.gz.md5

# hg38
if [ "$server" == "hpf" ]; then
    mv ${DIR_HG38}/clinvar.vcf.gz ${DIR_HG38}/clinvar_old/clinvar.${DATE}.vcf.gz
    rm ${DIR_HG38}/clinvar.vcf.gz.tbi ${DIR_HG38}/clinvar.vcf.gz.md5

    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz -O ${DIR_HG38}/clinvar.vcf.gz
    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi -O ${DIR_HG38}/clinvar.vcf.gz.tbi
    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.md5 -O ${DIR_HG38}/clinvar.vcf.gz.md5
    md5=`awk '{print $1}'  ${DIR_HG38}/clinvar.vcf.gz.md5`
    echo -e "${md5}\t${DIR_HG38}/clinvar.vcf.gz" > ${DIR_HG38}/clinvar.vcf.gz.md5
    md5sum -c ${DIR_HG38}/clinvar.vcf.gz.md5
fi
