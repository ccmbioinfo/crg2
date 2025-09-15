# crg2 - GRCm38 WES
Research pipeline for exploring variants in mouse exome (WES) sequencing data

<div align="center">
    <img src="/crg2logolarge.png" width="800px"</img> 
</div>

## Pipeline details
- Alignment to GRCm38/mm10 using BWA mem. 
- Duplicate marking with Picard. 
- Base quality score recalibration with GATK BaseRecalibrator. 
- Germline variant calling with GATK HaplotypeCaller.
- Somatic variant calling with GATK Mutect2, using MGP dbSNP142 variants as a germline resource. 
- Somatic variant QC with GATK PileupSummaries and GATK CalculateContamination.
- Somatic variant filtering with GATK FilterMutectCalls.
- Germline and somatic variant annotation using VEP. 

## Running the pipeline
1. Make a folder in a directory with sufficient space. Copy over the template files crg2/samples.tsv, crg2/units.tsv, crg2/config_hpf.yaml, crg2/dnaseq_slurm_hpf.sh, crg2/slurm_profile/slurm-config.yaml .
You may need to re-copy config_hpf.yaml and slurm-config.yaml if the files were recently updated in the repo from previous crg2 runs. Note that 'slurm-config.yaml' is for submitting each rule as cluster jobs, so ignore this if not running on cluster.
```
mkdir NA12878
cp crg2/samples.tsv crg2/units.tsv crg2/config_hpf.yaml crg2/dnaseq_slurm_hpf.sh crg2/slurm_profile/slurm-config.yaml NA12878
cd NA12878
```

2. Set up pipeline run 
  * Reconfigure 'samples.tsv', 'units.tsv' to reflect sample names and input files.
  * Modify 'config_hpf.yaml': 
    * change `project` to refer to the family ID (here, NA12878).

samples.tsv
```
sample
NA12878
```

units.tsv
```
sample	platform	fq1	fq2	bam	cram
NA12878	ILLUMINA			/hpf/largeprojects/ccmbio/GIAB_benchmark_datasets/hg19/WGS/NA12878/RMNISTHS_30xdownsample.bam
```

config_hpf.yaml
```
run:
  project: "NA12878"
  samples: samples.tsv
  units: units.tsv
  gatk: "gatk"
  pipeline: "wes"

...
```

3. Activate the conda environment with Snakemake 5.10.0

```
(base) [dennis.kao@qlogin5 crg2]$ conda activate snakemake
(snakemake) [dennis.kao@qlogin5 crg2]$ snakemake -v
5.10.0
```

4. Test that the pipeline will run by adding the flag "-n" to the command in dnaseq_slurm_hpf.sh and running it.


5. Submit pipeline job to the cluster

    Job scheduler: Slurm
      * Parallelized jobs on SickKids hpf:
      ```sbatch dnaseq_slurm_hpf.sh```
