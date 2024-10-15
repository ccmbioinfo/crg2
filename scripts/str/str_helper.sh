#!/bin/bash

x=`samtools idxstats $1 | grep "X"`;
y=`samtools idxstats $1 | grep "Y"`;
xcov=`echo $x | awk '{ printf("%0.5f", $3/$2); }'`;
ycov=`echo $y | awk '{ printf("%0.5f", $3/$2); }'`;

# xcov=$(echo "scale=4; $(samtools idxstats $1 | grep "X" | cut -f 3)/$(samtools idxstats $1 | grep "X" | cut -f 2)" | bc)
# ycov=$(echo "scale=4; $(samtools idxstats $1 | grep "Y" | cut -f 3)/$(samtools idxstats $1 | grep "Y" | cut -f 2)" | bc)
rat=$(echo "scale=4; ${xcov}/${ycov}" | bc)
if (( $(echo "$rat > 5.0" | bc -l) )); then
    sex=female
else
    sex=male
fi
echo "$sex"