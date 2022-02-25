#!/bin/bash 
#usage: ./copy_reports_to_minio.sh </path/to/family/results/dir> <NGS type: wgs or wes> <OPTIONAL: multiqc> <OPTIONAL: coverages>
 
family_dir=$1
NGS_type=$2
multiqc=$3
coverages=$4

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

if [ -z $NGS_type ]; then
        echo 'Must specify NGS type: wgs or wes'
        exit
fi

if [ "$NGS_type" != "wgs" ] && [ "$NGS_type" != "wes" ]; then
	echo 'NGS type must be wgs or wes, exiting'
	exit
fi

$mc mb minio/results-c4r/$family

if [ "$NGS_type" == wgs ]; then
	$mc cp ${family_dir}/report/coding/${family}/${family}.wes.regular*csv minio/results-c4r/$family
	$mc cp ${family_dir}/report/denovo/${family}/${family}.denovo*csv minio/results-c4r/$family
	$mc cp ${family_dir}/report/cnv/${family}*cnv.withSVoverlaps.tsv minio/results-c4r/$family
	$mc cp  ${family_dir}/report/sv/${family}.unfiltered.wgs*overlaps.tsv minio/results-c4r/$family
	$mc cp  ${family_dir}/report/sv/${family}.wgs*overlaps.tsv minio/results-c4r/$family
	$mc cp  ${family_dir}/report/str/${family}.EHDN.2021*xlsx minio/results-c4r/$family
	$mc cp  ${family_dir}/report/str/${family}.EH-v1.1.2021*xlsx minio/results-c4r/$family
	
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
else
	$mc cp ${family_dir}/${family}.clinical.wes.regular*csv minio/results-c4r/$family
	$mc cp ${family_dir}/${family}.clinical.wes.synonymous*csv minio/results-c4r/$family
	$mc cp ${family_dir}/${family}.wes.regular*csv minio/results-c4r/$family
	$mc cp ${family_dir}/${family}.wes.synonymous*csv minio/results-c4r/$family
fi
	
if [ "$multiqc" == multiqc ]; then
	$mc cp ${family_dir}/multiqc/multiqc_report.html minio/results-c4r/$family
fi

if [[ "$coverages" == coverages ]]; then
   if [ ! -d ${family_dir}/${family}_*_coverages ]; then
     echo 'No coverages directory found. Please make sure that your directory is named in the format of family_sample_coverages'
     exit
   fi
   for direc in ${family_dir}/${family}_*_coverages
     do 
     sample_cov=$(basename $direc)
     $mc mb minio/results-c4r/$family/$sample_cov
     $mc cp ${direc}/${family}*.dedup.bam.coverage_stats.csv minio/results-c4r/$family/$sample_cov
     $mc cp ${direc}/${family}*.dedup.bam.less20x_coverage.csv minio/results-c4r/$family/$sample_cov
     $mc cp ${direc}/${family}*.dedup.bam.less20x.stats.csv minio/results-c4r/$family/$sample_cov
     $mc cp ${direc}/${family}*.dedup.bam.median minio/results-c4r/$family/$sample_cov
     $mc cp ${direc}/${family}*.dedup.bam.per_exon.csv minio/results-c4r/$family/$sample_cov
     $mc cp ${direc}/${family}*.dedup.bam.per_exon.distribution.csv minio/results-c4r/$family/$sample_cov
     done
fi
