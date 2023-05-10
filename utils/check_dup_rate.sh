#!/bin/bash
family=$1
multiqc=`echo $family/qc/multiqc/multiqc_data/multiqc_picard_dups.txt`
<<<<<<< HEAD
dups=(`cat $multiqc | awk -F '\t' '{print $10}'`)
samples=(`cat $multiqc | awk -F '\t' '{print $1}'`)
num_of_samples=$((`echo ${#samples[@]}`-1))

for i in `seq 1 $num_of_samples`
do
	dup_rate=`awk "BEGIN {print ${dups[i]}*100}" | cut -d '.' -f1`
	sample=`echo ${samples[i]}`
	if (( $dup_rate >= "20")); then
		echo "$sample in $family has greater than 20% duplication rate, run coverage metrics: ${dups[i]}"
	else
		echo "$sample in $family duplication rate OK: ${dups[i]}"
=======
dups=`cat $multiqc | awk -F '\t' '{print $10}'`
for d in $dups
do
	if [ "$d" != "PERCENT_DUPLICATION" ]; then
		dups=`awk "BEGIN {print ${d}*100}" | cut -d '.' -f1`
		if (( $dups >= "20")); then
			echo "A sample in ${family} has greater than 20% duplication rate, run coverage metrics: $d"
		fi
>>>>>>> useful scripts for uploading/moving data and making bam slices
	fi
done

