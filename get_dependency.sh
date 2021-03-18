#!/bin/bash

#usage: sh get_dependency.sh
#generate this once; write a rule to copy the output file
#output=dependencies_<DATE>; run this from crg2 repo path


d=`date +%m-%d-%Y`;
out="programs.${d}.txt";
for i in `find -name "environment.yaml"`; do grep "=="  $i | awk '{ split($3,a,"=="); print $2,a[2];}' >> $out; done;
cat $out | sort -k1,1 -k2,2 | uniq >> temp && mv temp $out