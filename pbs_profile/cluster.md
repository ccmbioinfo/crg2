# Cluster integration

Snakemake can be made to submit each rule as a separate job to the PBS cluster by using `--profile__` or `_--cluster_` (deprecated). This is advantageous specifically on rules that can be run in parallel without waiting on inputs. In _crg2_, we have used the `--profile` option to make it more generic and manageable for future extension. 

https://snakemake.readthedocs.io/en/v5.1.4/executable.html#profiles
https://snakemake.readthedocs.io/en/v5.10.0/executing/cluster-cloud.html

```
pbs_profile is the profile folder with all the relevant files and settings for PBS job submission.
    |- config.yaml: contains common cluster parameters
    |- pbs_config.yaml: contains default and rule specific cluster settings like memory, nodes, walltime etc.,
    |- pbs_submit.py: script to create the 'qsub' commands for each job
    |- pbs_status.py: script to check/interpret job status
    |- pbs_jobscript.sh: template job script
    |- key_mapping.yaml: dictionary of various parameters available in 'qsub'
```
The above files for the profile were adapted from https://github.com/Snakemake-Profiles/generic and https://github.com/Snakemake-Profiles/pbs-torque

Copy the `pbs_config.yaml` file to the directory you wish to start Snakemake (this is expected to be in current path, just like the main crg2/config.yaml). Run Snakemake in dry-run mode by passing the absolute location of the profile folder as below,

```bash
snakemake --use-conda -s ~/crg2/Snakefile --conda-prefix /hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake --profile ~/crg2/pbs_profile -nr
```

If the profile directory is located in either of these two places `$HOME/.config/snakemake` and `/etc/xdg/snakemake`, then you can just pass the name of the profile folder. 

To submit snakemake to cluster, use `dnaseq_cluster.pbs` and change the value of variables inside as appropriate to your settings.


## Development notes

A. I have created group named 'gatkcall' encompassing the following rules, as this will group and submit single job for each contig without having to wait
    1. call_variants
    2. combine_calls
    3. genotype_variants

B. Setting `immediate_submit=true` in 'config.yaml' along with passing `--depend {dependencies}` to 'pbs_submit.py' will make Snakemake to submit all rules with job-chaining as per the DAG (`qsub -W "depend:afterok<jobids>"`) and exit. The problem with this approach is that, 
    1. when a job depends on multiple other jobs where all of them do not complete at the same time, some of the jobids are removed from the queue and 'qsub' will exit with submission error.  
    2. The main Snakemake process exits immediately, so there is no process to track output file presence before consequent job submission, monitor jobs and resubmit if failure happens. 
Therefore, setting 'immediate_submission=false' and not passing '--depend' will work best for our scenario. 

C.`pbs_submit.py` script prioritizes (mem and thread) options set in `pbs_config.yaml` over resources set inside ".smk" files. This would only affect the resources when submitting job, not the threads in actual program execution. May have to test for every case and rethink if this is not desired.  
