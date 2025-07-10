#!/bin/bash

# get OMIM files
set -e

cre_dir=/srv/shared/pipelines/cre
DATE=`date +%Y-%m-%d`
OMIM_key=$1 # e.g. lWLUD8lnRx-LhhI9832mXY

export PATH=:~/tools/:~/crg2/:~/bioscripts/:~/bioscripts/scripts/:/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/bin:~/cre/:/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/:$PATH

git clone https://github.com/ccmbioinfo/cre
cd cre
git checkout -b anno-update-${DATE}

# note: you must apply for OMIM download access once a year 
curl -o data/mim2gene.txt https://omim.org/static/omim/data/mim2gene.txt 
curl -o data/mimTitles.txt https://data.omim.org/downloads/${OMIM_key}/mimTitles.txt 
curl -o data/genemap2.txt https://data.omim.org/downloads/${OMIM_key}/genemap2.txt 
curl -o data/morbidmap.txt https://data.omim.org/downloads/${OMIM_key}/morbidmap.txt 

sh cre.update_omim.sh
previous_date=`grep 202 cre.vcf2db.R  | rev | cut -d '_' -f1 | rev | tr -d '.tsv")'`
sed -i "s/${previous_date}/${DATE}/g" cre.vcf2db.R

# get Orphanet file

sh cre.orphanet.sh
mv orphanet.txt data/orphanet.txt


git add data/hgnc_${DATE}.txt data/genemap2.txt data/mim2gene.txt data/mimTitles.txt data/OMIM_hgnc_join_omim_phenos_${DATE}.tsv data/OMIM_remaining_unmapped_omim_with_info_${DATE}.tsv data/orphanet.txt cre.vcf2db.R
git commit -m "Update OMIM and orphanet"
git push origin anno-update-${DATE}