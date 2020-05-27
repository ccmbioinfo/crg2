import pandas as pd


def get_contigs():
    return pd.read_table("/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa.fai",
                         header=None, usecols=[0], squeeze=True, dtype=str)


table = get_contigs()

print(table)
