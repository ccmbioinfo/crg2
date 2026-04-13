import pandas as pd

metadata = pd.read_table(snakemake.input.metadata, dtype=str)
metadata = metadata.replace(r"^\s*$", pd.NA, regex=True)
samples = pd.read_table(snakemake.input.samples, dtype=str)

required_cols = {"Sample", "Genotype", "Stage"}
missing = required_cols - set(metadata.columns)
if missing:
    raise ValueError(
        "Metadata file is missing required columns: {}".format(sorted(missing))
    )

if metadata["Sample"].duplicated().any():
    duplicated = sorted(metadata.loc[metadata["Sample"].duplicated(), "Sample"].tolist())
    raise ValueError("Duplicate Sample entries in metadata.tsv: {}".format(duplicated))

sample_set = set(samples["sample"])
metadata_set = set(metadata["Sample"])

missing_in_metadata = sorted(sample_set - metadata_set)
extra_in_metadata = sorted(metadata_set - sample_set)

if missing_in_metadata:
    raise ValueError(
        "Samples present in samples.tsv but missing from metadata.tsv: {}".format(missing_in_metadata)
    )

if extra_in_metadata:
    raise ValueError(
        "Samples present in metadata.tsv but not in samples.tsv: {}".format(extra_in_metadata)
    )

project_metadata = metadata[metadata["Sample"].isin(sample_set)].copy()
project_metadata.to_csv(snakemake.output.validated, sep="\t", index=False)
