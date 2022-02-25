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
echo $region

for bam in *.bam;
do
	if [ ! -d bam_slice ]; then
		mkdir bam_slice
	fi
	name=`echo $bam | cut -d'.' -f1`
	samtools view -b $bam "$region" > bam_slice/${name}_${region}.bam
	samtools index bam_slice/${name}_${region}.bam
done


