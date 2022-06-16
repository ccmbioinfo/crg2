#!/bin/bash

# reports holes in coverage: how many PB are below 20X per exon

echo "target,gene,pos,min_coverage,bases_below20" 

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
    if (prev_target != target && prev_target !=0){
	    print prev_target","prev_gene","prev_pos","prev_cov","prev_bases_below20;
	    prev_target=0;
    	    prev_pos=0;
    	    prev_cov=0;
    	    prev_gene=0;
    	    prev_bases_below20=0;
    }

    if ($4 < 20){
	prev_bases_below20+=1;
	if (prev_target !=0) {
	    if ($4 < prev_cov){
		prev_cov = $4
		prev_pos = $3
	    }
	}else{
	    prev_target=target;
    	    prev_pos=$3;
	    prev_cov=$4;
	    prev_gene=$2;
	}
    }
}
END{
    if (prev_target !=0){
	print prev_target","prev_gene","prev_pos","prev_cov","prev_bases_below20;
    }
}'
