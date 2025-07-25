#!/bin/bash 
#usage: ./copy_reports_to_minio.sh </path/to/family/results/dir> <NGS type: wgs or wes> <OPTIONAL: multiqc> <OPTIONAL: coverages>
 
family_dir=$1
NGS_type=$2
genome_type=$3
multiqc=$3
coverages=$4

mc=/hpf/largeprojects/ccmbio/ccmmarvin_shared/data_transfers/mc
if [ -z $family_dir ]; then
        echo 'Must specify results path, exiting'
        exit
fi

family=`echo $family_dir | cut -d '/' -f1`

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

if [ -z $genome_type ]; then
        echo 'Must specify genome type: TCAG or GSO'
        exit
fi

if [ "$genome_type" != "TCAG" ] && [ "$genome_type" != "GSO" ]; then
        echo 'NGS type must be TCAG or GSO, exiting'
        exit
fi

mkdir /srv/minio/results-c4r/$family

if [ "$NGS_type" == wgs ]; then
	cp ${family_dir}/report/coding/${family}/${family}.wes.regular*csv /srv/minio/results-c4r/$family
	cp ${family_dir}/report/denovo/${family}/${family}.denovo*csv /srv/minio/results-c4r/$family
	cp  ${family_dir}/report/sv/${family}.BND.wgs.sv*.tsv /srv/minio/results-c4r/$family
	cp  ${family_dir}/report/str/${family}.EHDN.202*xlsx /srv/minio/results-c4r/$family
	cp  ${family_dir}/report/str/${family}.EH-v1.1.2.202*xlsx /srv/minio/results-c4r/$family
	cp  ${family_dir}/report/mitochondrial/${family}.mitochondrial.report.202*csv /srv/minio/results-c4r/$family

	if [ "$genome_type" == TCAG ]; then
        	cp ${family_dir}/report/cnv/${family}*cnv.withSVoverlaps.tsv /srv/minio/results-c4r/$family
        	cp  ${family_dir}/report/sv/${family}.unfiltered.wgs*overlaps.tsv /srv/minio/results-c4r/$family
        	cp  ${family_dir}/report/sv/${family}.wgs*overlaps.tsv /srv/minio/results-c4r/$family
	else
		# no CNV reports
                cp  ${family_dir}/report/sv/${family}.unfiltered.wgs*tsv /srv/minio/results-c4r/$family
                cp  ${family_dir}/report/sv/${family}.wgs*tsv /srv/minio/results-c4r/$family
	fi
	
	if [ -f ${family_dir}/report/panel/${family}/${family}.wgs.panel*csv ]; then	
		cp ${family_dir}/report/panel/${family}/${family}.wgs.panel*csv /srv/minio/results-c4r/$family
	else
		panel=${family_dir}/report/panel/${family}/${family}.wgs*csv
		panel_date=`basename $panel | cut -d'.' -f3`
		mv $panel ${family_dir}/report/panel/${family}/${family}.wgs.panel.${panel_date}.csv
		cp ${family_dir}/report/panel/${family}/${family}.wgs.panel.${panel_date}.csv /srv/minio/results-c4r/$family
	fi 

	if [ -f ${family_dir}/report/panel-flank/${family}/${family}.wgs.panel-flank*csv ]; then
		cp ${family_dir}/report/panel-flank/${family}/${family}.wgs.panel-flank*csv /srv/minio/results-c4r/$family
	else
		panel_flank=${family_dir}/report/panel-flank/${family}/${family}.wgs*csv
		panel_flank_date=`basename $panel_flank | cut -d'.' -f3`
		mv $panel_flank ${family_dir}/report/panel-flank/${family}/${family}.wgs.panel-flank.${panel_flank_date}.csv
		cp ${family_dir}/report/panel-flank/${family}/${family}.wgs.panel-flank.${panel_flank_date}.csv /srv/minio/results-c4r/$family
	fi 
else
	cp ${family_dir}/report/coding/${family}/${family}.clinical.wes.regular*csv /srv/minio/results-c4r/$family
	cp ${family_dir}/report/coding/${family}/${family}.clinical.wes.synonymous*csv /srv/minio/results-c4r/$family
	cp ${family_dir}/report/coding/${family}/${family}.wes.regular*csv /srv/minio/results-c4r/$family
	cp ${family_dir}/report/coding/${family}/${family}.wes.synonymous*csv /srv/minio/results-c4r/$family
fi
	
if [ "$multiqc" == multiqc ]; then
	cp ${family_dir}/qc/multiqc/multiqc.html /srv/minio/results-c4r/$family
fi

if [[ "$coverages" == coverages ]]; then
   if [ ! -d ${family_dir}/${family}_*_coverages ]; then
     echo 'No coverages directory found. Please make sure that your directory is named in the format of family_sample_coverages'
     exit
   fi
   for direc in ${family_dir}/${family}_*_coverages
     do 
     sample_cov=$(basename $direc)
     mkdir /srv/minio/results-c4r/$family/$sample_cov
     cp ${direc}/${family}*.dedup.bam.coverage_stats.csv /srv/minio/results-c4r/$family/$sample_cov
     cp ${direc}/${family}*.dedup.bam.less20x_coverage.csv /srv/minio/results-c4r/$family/$sample_cov
     cp ${direc}/${family}*.dedup.bam.less20x.stats.csv /srv/minio/results-c4r/$family/$sample_cov
     cp ${direc}/${family}*.dedup.bam.median /srv/minio/results-c4r/$family/$sample_cov
     cp ${direc}/${family}*.dedup.bam.per_exon.csv /srv/minio/results-c4r/$family/$sample_cov
     cp ${direc}/${family}*.dedup.bam.per_exon.distribution.csv /srv/minio/results-c4r/$family/$sample_cov
     done
fi
