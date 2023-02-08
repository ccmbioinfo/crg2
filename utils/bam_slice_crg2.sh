#!/bin/bash
# usage: sh bam_slice.sh <region> <flank>
# example: sh bam_slice.sh  7:112725001-112728000 1000
# IMPORTANT: chromosome notation must match bam header (i.e. '7' for crg2 bams, 'chr7' for PacBio bams)
# run within <family>/recal 

region=$1
flank=$2
coord=`echo $region | cut -d':' -f2`
chr=`echo $region | cut -d':' -f1`
start_coord=`echo $coord | cut -d'-' -f1`
end_coord=`echo $coord | cut -d'-' -f2`
start_minus=`expr $start_coord - $flank`
end_minus=`expr $end_coord + $flank`
region=`echo ${chr}:${start_minus}-${end_minus}`
ref=/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa

echo $region

for cram in *.cram;
do
	if [ ! -d bam_slice ]; then
		mkdir bam_slice
	fi
	name=`echo $cram | cut -d'.' -f1`
	samtools view -b $cram "$region" -T $ref > bam_slice/${name}_${chr}-${start_minus}-${end_minus}.bam
	samtools index bam_slice/${name}_${chr}-${start_minus}-${end_minus}.bam
done

