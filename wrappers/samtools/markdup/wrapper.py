from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    family = snakemake.wildcards.family
    sample = snakemake.wildcards.sample

    ref_cache_cmd = "export REF_CACHE={snakemake.params.ref_cache}; export REF_PATH={snakemake.params.ref_cache}; "
    name_sort_cmd = "samtools sort -n -o dups_removed/{family}_{sample}_namesort.cram {snakemake.input.cram}; "
    fixmate_cmd = "samtools fixmate -m dups_removed/{family}_{sample}_namesort.cram dups_removed/{family}_{sample}_fixmate.cram; "
    pos_sort_cmd = "samtools sort -o  dups_removed/{family}_{sample}_fixmate_positionsort.cram dups_removed/{family}_{sample}_fixmate.cram; "
    rm_dup_cmd = "samtools markdup -r dups_removed/{family}_{sample}_fixmate_positionsort.cram dups_removed/{family}_{sample}.cram; "
    rm_intermediates_cmd = "rm dups_removed/{family}_{sample}_namesort.cram dups_removed/{family}_{sample}_fixmate_positionsort.cram dups_removed/{family}_{sample}_fixmate.cram"

    shell(
    "("
    + ref_cache_cmd
    + name_sort_cmd
    + fixmate_cmd
    + pos_sort_cmd
    + rm_dup_cmd
    + rm_intermediates_cmd
    + ")"
    )
