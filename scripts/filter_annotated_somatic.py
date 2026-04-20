import pysam

family = snakemake.wildcards.family
MIN_DP = snakemake.params.min_dp
MIN_AF = snakemake.params.min_af
MIN_ALT_AC = snakemake.params.min_alt_ac

somatic_samples = snakemake.params.somatic_samples
somatic_vcf_samples = [f"{family}_{sample}" for sample in somatic_samples]

vcf_in = pysam.VariantFile(snakemake.input.vcf)
vcf_out = pysam.VariantFile(snakemake.output.vcf, "wz", header=vcf_in.header)

for record in vcf_in:
    keep = False

    for sample in somatic_vcf_samples:
        sample_data = record.samples[sample]

        dp = sample_data.get("DP")
        dp_ok = dp is not None and dp >= MIN_DP

        af = sample_data.get("AF")
        af_ok = False
        if af is not None:
            if isinstance(af, tuple):
                af_ok = any(x is not None and x >= MIN_AF for x in af)
            else:
                af_ok = af >= MIN_AF

        ad = sample_data.get("AD")
        alt_ac_ok = False
        if ad is not None and len(ad) > 1:
            alt_ac_ok = any(x is not None and x >= MIN_ALT_AC for x in ad[1:])

        if dp_ok and af_ok and alt_ac_ok:
            keep = True
            break

    if keep:
        vcf_out.write(record)

vcf_in.close()
vcf_out.close()
