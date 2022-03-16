
#!/bin/bash

read -r -a l <<< $@;
a=($(echo ${l[4]} | awk '{ split($1,x,""); print x[1],x[2],x[3],x[4],x[5];}')); 
callers=(); 
sum=$(IFS=+; echo "$((${a[*]}))")
if [[ ${a[0]} == 1 ]]; then callers+=("gatk-haplotype"); fi;  
if [[ ${a[1]} == 1 ]]; then callers+=("samtools"); fi; 
if [[ ${a[2]} == 1 ]]; then callers+=("freebayes"); fi; 
if [[ ${a[3]} == 1 ]]; then callers+=("platypus"); fi;
if [[ ${a[4]} == 1 ]]; then callers+=("freebayes-mosaic"); fi;  
callers=$(IFS=","; echo "${callers[*]};");  
echo -e "${l[0]}\t${l[1]}\t${l[2]}\t${l[3]}\t${callers}\t${sum}"; 

