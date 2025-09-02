#!/bin/bash
# Update OMIM and Orphanet annotations for cre, the cre-hg38 branch, and crg2-pacbio cre
#
# This script will:
# 1. Update OMIM and Orphanet annotations in both repos
# 2. Create new branches and push changes to remote repos
# 3. After running, you'll need to create and merge PRs in each repo
#
# Usage:
#   sh update_annotations.sh <OMIM_key>
#
# Arguments:
#   OMIM_key    Access key from OMIM (e.g. lWLUD8lnRx-LhhI9832mXY)
#               You must apply for OMIM download access annually.
#               The key can be found in the download links they email you.

# get OMIM files
set -e

DATE=`date +%Y-%m-%d`
OMIM_key=$1 # e.g. lWLUD8lnRx-LhhI9832mXY

export PATH=:~/tools/:~/crg2/:~/bioscripts/:~/bioscripts/scripts/:/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/bin:~/cre/:/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/:$PATH # for Rscript

# clone fresh instance of cre repo
git clone https://github.com/ccmbioinfo/cre
cd cre
git checkout -b anno-update-${DATE}

# download OMIM files
curl -o data/mim2gene.txt https://omim.org/static/omim/data/mim2gene.txt 
curl -o data/mimTitles.txt https://data.omim.org/downloads/${OMIM_key}/mimTitles.txt 
curl -o data/genemap2.txt https://data.omim.org/downloads/${OMIM_key}/genemap2.txt 
curl -o data/morbidmap.txt https://data.omim.org/downloads/${OMIM_key}/morbidmap.txt 

sh cre.update_omim.sh `pwd`

# update cre.vcf2db.R to use new OMIM files
previous_date=`grep 202 cre.vcf2db.R  | rev | cut -d '_' -f1 | rev | tr -d '.tsv")'`
sed -i "s/${previous_date}/${DATE}/g" cre.vcf2db.R

# get Orphanet file

sh cre.orphanet.sh
mv orphanet.txt data/orphanet.txt

# now add all changes to the cre branch and commit
git add data/hgnc_${DATE}.txt data/genemap2.txt data/mim2gene.txt data/mimTitles.txt data/morbidmap.txt data/OMIM_hgnc_join_omim_phenos_${DATE}.tsv data/OMIM_remaining_unmapped_omim_with_info_${DATE}.tsv data/orphanet.txt cre.vcf2db.R
git commit -m "Update OMIM and orphanet"
#git push origin anno-update-${DATE}

# update same files for cre hg38 branch
cd ../
mv cre cre-hg19
git clone https://github.com/ccmbioinfo/cre
cd cre
git checkout hg38
git checkout -b anno-update-${DATE} hg38

# copy OMIM files
cp ../cre-hg19/data/OMIM_hgnc_join_omim_phenos_${DATE}.tsv ../cre-hg19/data/OMIM_remaining_unmapped_omim_with_info_${DATE}.tsv data
cp ../cre-hg19/data/genemap2.txt ../cre-hg19/data/mim2gene.txt ../cre-hg19/data/mimTitles.txt ../cre-hg19/data/hgnc_${DATE}.txt  data
cp ../cre-hg19/data/orphanet.txt data 

# update cre.vcf2db.R to use new OMIM files
previous_date=`grep 202 cre.vcf2db.R  | rev | cut -d '_' -f1 | rev | tr -d '.tsv")'`
sed -i "s/${previous_date}/${DATE}/g" cre.vcf2db.R

# now add all changes to the cre branch and commit
git add data/hgnc_${DATE}.txt data/genemap2.txt data/mim2gene.txt data/mimTitles.txt data/OMIM_hgnc_join_omim_phenos_${DATE}.tsv data/OMIM_remaining_unmapped_omim_with_info_${DATE}.tsv data/orphanet.txt cre.vcf2db.R
git commit -m "Update OMIM and orphanet"
git push origin anno-update-${DATE}


# update same files for crg2-pacbio cre
cd ../
git clone https://github.com/ccmbioinfo/crg2-pacbio
cd crg2-pacbio
git checkout -b anno-update-${DATE}

cp ../cre-hg19/data/OMIM_hgnc_join_omim_phenos_${DATE}.tsv  scripts/cre/data
cp ../cre-hg19/data/genemap2.txt ../cre-hg19/data/mim2gene.txt ../cre-hg19/data/mimTitles.txt ../cre-hg19/data/morbidmap.txt  scripts/cre/data
cp ../cre-hg19/data/orphanet.txt scripts/cre/data 

# update cre.vcf2db.R
previous_date=`grep 202 scripts/cre/cre.vcf2db.R  | rev | cut -d '_' -f1 | rev | tr -d '.tsv")'`
sed -i "s/${previous_date}/${DATE}/g" scripts/cre/cre.vcf2db.R

# push changes to crg2-pacbio repo
git add scripts/cre/data/OMIM_hgnc_join_omim_phenos_${DATE}.tsv scripts/cre/data/orphanet.txt scripts/cre/cre.vcf2db.R scripts/cre/data/genemap2.txt scripts/cre/data/mim2gene.txt scripts/cre/data/mimTitles.txt scripts/cre/data/morbidmap.txt
git commit -m "Update OMIM and orphanet"
git push origin anno-update-${DATE}

