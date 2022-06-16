#!/bin/env python

# $1 = file.bam.coverage = product of ~/bioscripts/bam.coverage.sh with -d on
# $2 = No of the field with coverage, 0-based
# $3 = character separator in the input file

import sys
import numpy

cov={}

levels = [0,1,10,20,30,50,100,200]
for i in levels:
    cov[i]=0

total_bases=0

coverage_field_num = int(sys.argv[2])
sep = sys.argv[3]

with open(sys.argv[1],'rb') as f_coverage:
    for line in f_coverage:
	total_bases += 1
	buf = line.strip()
	fields = buf.split(sep)
	coverage = int(fields[coverage_field_num])

	for j in levels:
	    if (coverage >= j):
		cov[j] += 1

f_coverage.close()

print("Coverage,bases at,%")

for i in levels:
    print(str(i)+','+str(cov[i])+','+str(1.0*cov[i]/total_bases))

