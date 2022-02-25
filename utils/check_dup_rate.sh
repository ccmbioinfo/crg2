#!/bin/bash
family=$1
multiqc=`echo $family/qc/multiqc/multiqc_data/multiqc_picard_dups.txt`
dups=`cat $multiqc | awk -F '\t' '{print $10}'`
for d in $dups
do
	if [ "$d" != "PERCENT_DUPLICATION" ]; then
		dups=`awk "BEGIN {print ${d}*100}" | cut -d '.' -f1`
		if (( $dups >= "20")); then
			echo "A sample in ${family} has greater than 20% duplication rate, run coverage metrics: $d"
		fi
	fi
done

