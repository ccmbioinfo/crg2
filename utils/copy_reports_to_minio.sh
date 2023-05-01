#!/bin/bash 
#usage: ./copy_reports_to_minio.sh </path/to/family/results/dir> <NGS type: wgs or wes> <OPTIONAL: multiqc> <OPTIONAL: coverages>
 
family_dir=$1
genome_type=$2

mc=/hpf/largeprojects/ccmbio/ccmmarvin_shared/data_transfers/mc
if [ -z $family_dir ]; then
        echo 'Must specify results path, exiting'
        exit
fi

family=`basename $family_dir`

if [ ! -d $family_dir ]; then
	echo 'Directory does not exist, exiting'
	exit
fi

if [ -z $genome_type ]; then
        echo 'Must specify genome type: TCAG or GSO'
        exit
fi


$mc mb minio/results-c4r/$family

$mc cp ${family_dir}/report/coding/${family}/${family}.wes.regular*csv minio/results-c4r/$family
$mc cp ${family_dir}/report/denovo/${family}/${family}.denovo*csv minio/results-c4r/$family
$mc cp ${family_dir}/report/cnv/${family}*cnv.withSVoverlaps.tsv minio/results-c4r/$family
$mc cp  ${family_dir}/report/sv/${family}.unfiltered.wgs*overlaps.tsv minio/results-c4r/$family
$mc cp  ${family_dir}/report/sv/${family}.wgs*overlaps.tsv minio/results-c4r/$family
$mc cp  ${family_dir}/report/str/${family}.EHDN.202*xlsx minio/results-c4r/$family
$mc cp  ${family_dir}/report/str/${family}.EH-v1.1.202*xlsx minio/results-c4r/$family
	
if [ $genome_type == TCAG ]; then
	

if [ -f ${family_dir}/report/panel/${family}/${family}.wgs.panel*csv ]; then	
	$mc cp ${family_dir}/report/panel/${family}/${family}.wgs.panel*csv minio/results-c4r/$family
else
	panel=${family_dir}/report/panel/${family}/${family}.wgs*csv
	panel_date=`basename $panel | cut -d'.' -f3`
	mv $panel ${family_dir}/report/panel/${family}/${family}.wgs.panel.${panel_date}.csv
	$mc cp ${family_dir}/report/panel/${family}/${family}.wgs.panel.${panel_date}.csv minio/results-c4r/$family
fi 

if [ -f ${family_dir}/report/panel-flank/${family}/${family}.wgs.panel-flank*csv ]; then
	$mc cp ${family_dir}/report/panel-flank/${family}/${family}.wgs.panel-flank*csv minio/results-c4r/$family
else
	panel_flank=${family_dir}/report/panel-flank/${family}/${family}.wgs*csv
	panel_flank_date=`basename $panel_flank | cut -d'.' -f3`
	mv $panel_flank ${family_dir}/report/panel-flank/${family}/${family}.wgs.panel-flank.${panel_flank_date}.csv
	$mc cp ${family_dir}/report/panel-flank/${family}/${family}.wgs.panel-flank.${panel_flank_date}.csv minio/results-c4r/$family
fi 
	

