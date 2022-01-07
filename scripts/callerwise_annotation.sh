#!/bin/bash 

#usage: sh callerwise_annotation.sh <caller>[gatk-haplotype|samtools|platypus|freebayes] <sites.txt> <outdir> <out_vcf>

caller=$1;
sites=$2;
outdir=$3;
out_vcf=$4;

echo "extracting variants for $caller"

if [ ${caller} == "platypus" ]; then
grep "${caller}" ${sites} | egrep -v "gatk-haplotype|samtools|freebayes" > ${outdir}/${caller}.annot.txt
file=${outdir}/0003.vcf.gz

elif [ ${caller} == "freebayes" ]; then
file=${outdir}/0002.vcf.gz
grep "${caller}" ${sites} | egrep -v "gatk-haplotype|samtools" > ${outdir}/${caller}.annot.txt

elif [ ${caller} == "samtools" ]; then
grep ${caller} ${sites} | egrep -v "gatk-haplotype" > ${outdir}/${caller}.annot.txt
file=${outdir}/0001.vcf.gz

elif [ ${caller} == "gatk-haplotype" ]; then
grep ${caller} ${sites} > ${outdir}/${caller}.annot.txt
file=${outdir}/0000.vcf.gz

else
echo "${caller} caller not recognized";

fi;
bgzip -c ${outdir}/${caller}.annot.txt > ${outdir}/${caller}.annot.txt.gz
tabix -s1 -b2 -e2 ${outdir}/${caller}.annot.txt.gz

echo "annotating caller info to $file -> ${out_vcf}"

bcftools annotate -a ${outdir}/${caller}.annot.txt.gz -h ${outdir}/hdr.txt -c CHROM,POS,REF,ALT,INFO/CALLERS,INFO/NUMCALLS --threads 1  -o ${out_vcf} -O z $file
tabix ${out_vcf}

