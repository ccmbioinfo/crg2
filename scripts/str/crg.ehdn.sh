#PBS -N EHDN
#PBS -l vmem=80g,mem=80g,walltime=24:00:00
#PBS -joe
#PBS -d .

if [ -z $EHDN ]
then
    exit;
fi

if [ -z ${g1k_manifest} ]
then
    exit;
fi

if [ -z $ref ]
then
    exit;
fi

if [ -z $script_dir ]
then
    exit;
fi

scripts="${script_dir}";



if [ -z $family ]; then family=$1; fi;
if [ -z $pipeline ]; then pipeline=$2; fi #pipeline is crg or crg2

if [ $pipeline == "crg" ]; then 
    bcbio="${family}/bcbio-align/${family}/final/${family}_*";
    if [ -z "`ls ${bcbio}/${family}_*-ready.bam`" ]; then 
        if [ -z "`ls ${family}/${family}*-ready.bam`" ]; then
            echo "BAM files for $family not found inside bcbio-align/ and $family/. Exiting!";
            exit
        else
            dir="$family/${family}*-ready.bam";
        fi;
    else
        dir="$bcbio/${family}*-ready.bam";
    fi;
    outdir="${family}/str/expansion_hunter_denovo";
else
    dir="decoy_rm/*.bam"
    outdir="str/EHDN";
fi;

manifest="${outdir}/${family}_manifest.txt";
echo "Running for pipeline: ${pipeline} with directories => ${outdir}, ${dir}, $family";

for i in `ls $dir`; do

    if [ ! -d "$outdir" ]; then
        mkdir -p $outdir;
    fi;
    if [ $pipeline == "crg" ]; then 
        prefix=`basename $i -ready.bam`;
    else
        prefix=`basename $i .bam | cut -d "." -f1`; 
    fi
    reads=`readlink -f $i`;
    ehdn_prefix="${outdir}/${prefix}";
    echo -e "STEP1: ExpansionHunterDenovo-v0.7.0 --reference $ref --reads $reads --output-prefix $ehdn_prefix --min-anchor-mapq 50 --max-irr-mapq 40 \n"
    $EHDN --reference $ref --reads $reads --output-prefix $ehdn_prefix --min-anchor-mapq 50 --max-irr-mapq 40 &
    echo -e "${prefix}\tcase\t${ehdn_prefix}.json" >> $manifest;

done;

wait

cat $manifest ${g1k_manifest} > temp && mv temp $manifest;

##EHDN report per family
echo "STEP2: generating multi-sample EHDN STR profile and report named ${outdir}/${family}_EHDN_str.tsv"
#module load python/3.7.7

countsfile="${outdir}/${family}_counts.txt";
profile="${outdir}/${family}_EHDN_str.tsv";
python ${scripts}/combine_counts.py --manifest $manifest --combinedCounts $countsfile
python ${scripts}/compare_anchored_irrs.py --manifest $manifest --inputCounts $countsfile --outputRegions $profile --minCount 2 --testMethod normal
#module unload python/3.7.7

##add chr to chromosome column
echo "STEP2.1: add chr prefix to first column ${outdir}/${family}_EHDN_str.tsv"
awk -vFS="\t" -vOFS="\t"  '{ if(!($1~/^chr/)) {$1="chr"$1;} print $0; }' ${profile} > temp && mv temp ${profile}
