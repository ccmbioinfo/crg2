import pandas as pd

def parse_config_csv(value):
    if isinstance(value, str):
        return {item.strip() for item in value.split(",") if item.strip()}
    return set(value)

metadata = pd.read_table(snakemake.input.metadata, dtype=str)

baseline_genotype = str(snakemake.config["baseline"]["genotype"]).strip()
baseline_stages = parse_config_csv(snakemake.config["baseline"]["stages"])

baseline_mask = (
    (metadata["Genotype"] == baseline_genotype)
    & (metadata["Stage"].isin(baseline_stages))
)

baseline_samples = metadata.loc[baseline_mask, "Sample"].dropna().tolist()
nonbaseline_samples = metadata.loc[~baseline_mask, "Sample"].dropna().tolist()

with open(snakemake.output.baseline, "w") as handle:
    for sample in baseline_samples:
        handle.write(f"{sample}\n")

with open(snakemake.output.nonbaseline, "w") as handle:
    for sample in nonbaseline_samples:
        handle.write(f"{sample}\n")
