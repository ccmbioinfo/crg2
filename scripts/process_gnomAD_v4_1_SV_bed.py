import pandas as pd


gnomad_df = pd.read_csv("/hpf/largeprojects/ccmbio/ajain/crg2_hg38/gnomad_v4.1/gnomad_SV/gnomad.v4.1.sv.sites.bed", sep="\t", dtype="str").astype(str)
gnomad_df.columns = gnomad_df.columns.str.replace("#", "")
gnomad_df.columns = gnomad_df.columns.str.strip()
gnomad_df.rename({"chrom": "CHROM", "start": "START", "end": "END", "name": "NAME", "svtype": "SVTYPE", "samples": "SAMPLES" }, axis=1, inplace=True)
gnomad_df["CHROM"] = gnomad_df["CHROM"].str.replace("chr", "")
# in gnomAD v4.1 bed file, 'END' and 'SVTYPE' columns are duplicated
gnomad_df = gnomad_df.loc[:, ~gnomad_df.columns.duplicated()]
# select only columns of interest
gnomad_cols = [
    "CHROM",
    "START",
    "END",
    "NAME",
    "SVTYPE",
    "AN",
    "AC",
    "AF",
    "N_HOMREF",
    "N_HET",
    "N_HOMALT",
    "FREQ_HOMREF",
    "FREQ_HET",
    "FREQ_HOMALT",
    "GRPMAX_AF",
    "FILTER",
]
gnomad_df = gnomad_df[gnomad_cols]
gnomad_df.to_csv("/hpf/largeprojects/ccmbio/ajain/crg2_hg38/gnomad_v4.1/gnomad_SV/gnomad.v4.1.sv.sites.processed.bed", sep="\t", index=False )