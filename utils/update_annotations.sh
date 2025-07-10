#!/bin/bash

# get OMIM files
set -e

cre_dir=~/cre
OMIM_key=$1 # e.g. lWLUD8lnRx-LhhI9832mXY

export PATH=:~/tools/:~/crg2/:~/bioscripts/:~/bioscripts/scripts/:/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/bin:~/cre/:/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/:$PATH


# note: you must apply for OMIM download access once a year 
curl -o ${cre_dir}/data/mim2gene.txt https://omim.org/static/omim/data/mim2gene.txt 
curl -o ${cre_dir}/data/mimTitles.txt https://data.omim.org/downloads/${OMIM_key}/mimTitles.txt 
curl -o ${cre_dir}/data/genemap2.txt https://data.omim.org/downloads/${OMIM_key}/genemap2.txt 
curl -o ${cre_dir}/data/morbidmap.txt https://data.omim.org/downloads/${OMIM_key}/morbidmap.txt 

sh ~/cre/cre.update_omim.sh

# get Orphanet file

sh ~/cre/cre.orphanet.sh
mv orphanet.txt ${cre_dir}/data/orphanet.txt




