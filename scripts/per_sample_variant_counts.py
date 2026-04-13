import pysam
import pandas as pd

GERMLINE_MIN_DP = snakemake.params.germline_min_dp
GERMLINE_MIN_ALT_DEPTH = snakemake.params.germline_min_alt_depth

SOMATIC_MIN_DP = snakemake.params.somatic_min_dp
SOMATIC_MIN_AF = snakemake.params.somatic_min_af
SOMATIC_MIN_ALT_AC = snakemake.params.somatic_min_alt_ac

vcfs = {
    "germline_annotated": snakemake.input.germline_annotated,
    "germline_baseline_excluded": snakemake.input.germline_baseline_excluded,
    "germline_final": snakemake.input.germline_final,
    "somatic_annotated": snakemake.input.somatic_annotated,
    "somatic_baseline_excluded": snakemake.input.somatic_baseline_excluded,
    "somatic_final": snakemake.input.somatic_final,
}

rows = []


def has_nonref_gt(sample_data):
    gt = sample_data.get("GT")
    if gt is None:
        return False
    alleles = [a for a in gt if a is not None]
    if len(alleles) == 0:
        return False
    return any(a != 0 for a in alleles)


def get_max_alt_depth(sample_data):
    ad = sample_data.get("AD")
    if ad is None or len(ad) <= 1:
        return 0
    alt_depths = [x for x in ad[1:] if x is not None]
    return max(alt_depths) if alt_depths else 0


def passes_germline(sample_data):
    if not has_nonref_gt(sample_data):
        return False

    dp = sample_data.get("DP")
    if dp is None or dp < GERMLINE_MIN_DP:
        return False

    alt_depth = get_max_alt_depth(sample_data)
    if alt_depth < GERMLINE_MIN_ALT_DEPTH:
        return False

    return True


def passes_somatic(sample_data):
    if not has_nonref_gt(sample_data):
        return False

    dp = sample_data.get("DP")
    if dp is None or dp < SOMATIC_MIN_DP:
        return False

    alt_depth = get_max_alt_depth(sample_data)
    if alt_depth < SOMATIC_MIN_ALT_AC:
        return False

    af = sample_data.get("AF")
    if af is None:
        return False

    if isinstance(af, tuple):
        af_values = [x for x in af if x is not None]
        if not af_values or max(af_values) < SOMATIC_MIN_AF:
            return False
    else:
        if af < SOMATIC_MIN_AF:
            return False

    return True


for stage, path in vcfs.items():
    vcf = pysam.VariantFile(path)
    samples = list(vcf.header.samples)

    gt_counts = {sample: 0 for sample in samples}
    filtered_counts = {sample: 0 for sample in samples}

    is_germline = stage.startswith("germline_")

    for record in vcf:
        for sample in samples:
            sample_data = record.samples[sample]

            if has_nonref_gt(sample_data):
                gt_counts[sample] += 1

            if is_germline:
                if passes_germline(sample_data):
                    filtered_counts[sample] += 1
            else:
                if passes_somatic(sample_data):
                    filtered_counts[sample] += 1

    vcf.close()

    for sample in samples:
        rows.append(
            {
                "family": snakemake.wildcards.family,
                "sample": sample,
                "stage": stage,
                "n_variants_gt": gt_counts[sample],
                "n_variants_filtered": filtered_counts[sample],
            }
        )

pd.DataFrame(rows).to_csv(snakemake.output.summary, sep="\t", index=False)
