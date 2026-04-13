import pysam

family = snakemake.wildcards.family
MIN_DP = snakemake.params.min_dp
MIN_ALT_DEPTH = snakemake.params.min_alt_depth

nonbaseline_samples = snakemake.params.nonbaseline_samples
nonbaseline_vcf_samples = [f"{family}_{sample}" for sample in nonbaseline_samples]

vcf_in = pysam.VariantFile(snakemake.input.vcf)
vcf_out = pysam.VariantFile(snakemake.output.vcf, "wz", header=vcf_in.header)

for record in vcf_in:
    keep = False

    for sample in nonbaseline_vcf_samples:
        sample_data = record.samples[sample]

        dp = sample_data.get("DP")
        dp_ok = dp is not None and dp >= MIN_DP

        ad = sample_data.get("AD")
        alt_depth_ok = False
        if ad is not None and len(ad) > 1:
            alt_depth_ok = any(x is not None and x >= MIN_ALT_DEPTH for x in ad[1:])

        if dp_ok and alt_depth_ok:
            keep = True
            break

    if keep:
        vcf_out.write(record)

vcf_in.close()
vcf_out.close()
