import pysam
from collections import defaultdict

GERMLINE_MIN_DP = snakemake.params.germline_min_dp
GERMLINE_MIN_ALT_DEPTH = snakemake.params.germline_min_alt_depth


def has_nonref_gt(sample_data):
    gt = sample_data.get("GT")
    if gt is None:
        return False
    alleles = [a for a in gt if a is not None]
    if len(alleles) == 0:
        return False
    return any(a != 0 for a in alleles)


def max_alt_depth(sample_data):
    ad = sample_data.get("AD")
    if ad is None or len(ad) <= 1:
        return 0
    alt_depths = [x for x in ad[1:] if x is not None]
    return max(alt_depths) if alt_depths else 0


def passes_germline_threshold(sample_data):
    if not has_nonref_gt(sample_data):
        return False

    dp = sample_data.get("DP")
    if dp is None or dp < GERMLINE_MIN_DP:
        return False

    if max_alt_depth(sample_data) < GERMLINE_MIN_ALT_DEPTH:
        return False

    return True


def record_key(record):
    return (
        record.contig,
        record.pos,
        record.ref,
        tuple(record.alts or ())
    )


def zero_gt(gt):
    if gt is None:
        return None
    return tuple(0 if allele is not None else None for allele in gt)


def zero_ad(dp, n_alts, fallback_ref=0):
    ref_depth = dp if dp is not None else fallback_ref
    return tuple([ref_depth] + [0] * n_alts)


def zero_af(n_alts):
    return tuple([0.0] * n_alts)


germline_sites = defaultdict(set)

germline_vcf = pysam.VariantFile(snakemake.input.germline_vcf)
for record in germline_vcf:
    key = record_key(record)
    for sample in germline_vcf.header.samples:
        sample_data = record.samples[sample]
        if passes_germline_threshold(sample_data):
            germline_sites[sample].add(key)
germline_vcf.close()

somatic_vcf = pysam.VariantFile(snakemake.input.somatic_vcf)
out_vcf = pysam.VariantFile(snakemake.output.vcf, "wz", header=somatic_vcf.header)

for record in somatic_vcf:
    key = record_key(record)
    n_alts = len(record.alts or ())

    for sample in somatic_vcf.header.samples:
        if key not in germline_sites.get(sample, set()):
            continue

        sample_data = record.samples[sample]

        gt = sample_data.get("GT")
        dp = sample_data.get("DP")
        ad = sample_data.get("AD")
        af = sample_data.get("AF")

        if gt is not None:
            sample_data["GT"] = zero_gt(gt)

        if ad is not None:
            fallback_ref = ad[0] if len(ad) > 0 and ad[0] is not None else 0
            sample_data["AD"] = zero_ad(dp, n_alts, fallback_ref)

        if af is not None:
            sample_data["AF"] = zero_af(n_alts)

    keep = False
    for sample in somatic_vcf.header.samples:
        if has_nonref_gt(record.samples[sample]):
            keep = True
            break

    if keep:
        out_vcf.write(record)

somatic_vcf.close()
out_vcf.close()
