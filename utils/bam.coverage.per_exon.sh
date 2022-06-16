#!/bin/bash

# reports the distribution of lowest coverage in exons
# each exon is characterized by a bp with the lowest coverage
# distribution of theses values for 
# input is dcoverage

echo "target,gene,min_coverage"

cat $1 | awk -F "\t" '{print $1"-"$2"-"$3","$4","$2+$5","$6}' | \
awk -F ',' '
BEGIN{
    prev_target=0;
    prev_pos=0;
    prev_cov=0;
    prev_gene=0;
    prev_bases_below20=0;
}
{
    target=$1
    if (prev_target != target){
	    if (prev_target !=0){
		print prev_target","prev_gene","prev_cov;
	    }
	    prev_target=target;
    	    prev_cov=$4;
    	    prev_gene=$2;
    }
    if ($4 < prev_cov){
		prev_cov = $4
    }
}
END{
    print prev_target","prev_gene","prev_cov;
}'

