import pandas as pd
import pybedtools
import logging, os


def log_message(*message):
    """write message to logfile and stdout"""
    if message:
        for i in message:
            logging.info(i)
            print(i)


def query_alias(aliases, gene):
    """check if gene is in list of aliases"""
    # check if value is nan
    if aliases == aliases:
        aliases = [alias.replace(" ", "") for alias in aliases]
        if gene in aliases:
            return True
        else:
            return False
    else:
        return False


def main(family, hpo, ensembl, refseq, hgnc):
    logfile = "logs/hpo_to_panel/genes.log"
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.DEBUG,
        format="%(asctime)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M",
    )

    log_message("Loading gene coordinates and aliases")
    hpo = pd.read_csv(hpo, sep="\t", comment="#").set_index("Gene ID")
    ensembl = pd.read_csv(ensembl).astype(str).set_index("name")
    refseq = pd.read_csv(refseq).astype(str).set_index("name")
    hgnc = pd.read_csv(hgnc, sep="\t")
    genes = pd.concat([ensembl, refseq])

    # join hpo terms and ensembl+refseq genes on gene id (most gene ids in the hpo file are ENSG ids, but not all)
    log_message("Mapping HPO genes to gene coordinates")
    hpo_coord = genes.join(hpo, how="inner").reset_index()
    # add 'chr' for hg38 pipeline
    hpo_coord['chromosome'] = 'chr' + hpo_coord['chromosome'].astype(str)
    missing = [gene for gene in hpo.index if gene not in hpo_coord["name"].values]

    # check missing genes against HGNC previous symbols and aliases
    log_message("Mapping missing genes to HGNC aliases")
    hgnc["Alias symbols"] = hgnc["Alias symbols"].apply(
        lambda x: x.split(",") if x == x else x
    )
    hgnc["Previous symbols"] = hgnc["Previous symbols"].apply(
        lambda x: x.split(",") if x == x else x
    )
    alias_list = hgnc["Alias symbols"].dropna().tolist()
    previous_list = hgnc["Previous symbols"].dropna().tolist()

    alias_prev_symbols = alias_list + previous_list
    alias_prev_symbols = [
        alias.replace(" ", "") for alias in alias_prev_symbols for alias in alias
    ]

    found = []
    for gene in missing:
        if gene in alias_prev_symbols:
            try:
                symbol = hgnc[
                    hgnc["Alias symbols"].apply(lambda x: query_alias(x, gene))
                ]["Approved symbol"].values[0]
            except IndexError:
                symbol = hgnc[
                    hgnc["Previous symbols"].apply(lambda x: query_alias(x, gene))
                ]["Approved symbol"].values[0]
            # using approved symbol, check refseq genes again
            # if gene is in refseq, add coordinates to hpo_coord df
            try:
                coord = refseq.loc[symbol].copy()
                hpo_coord = pd.concat([hpo_coord, coord], axis=0)
                found.append(gene)
            except KeyError:
                continue

    missing = [gene for gene in missing if gene not in found]
    if len(missing) != 0:
        log_message(f"The following genes could not be found: {missing}")
        missing_df = pd.DataFrame(missing, columns=["Gene ID"])
        missing_df.to_csv("genes/missing_genes.txt", sep="\t", index=False)

    # subset columns so dataframe conforms to bed format
    hpo_coord = hpo_coord[["chromosome", "start", "end"]]
    hpo_coord = hpo_coord.dropna()
    hpo_coord.to_csv("genes/temp.bed", sep="\t", index=False, header=False)

    # create bedtools object for merging and sorting
    bed = pybedtools.BedTool("genes/temp.bed")
    bed = bed.sort().merge()
    bed_df = bed.to_dataframe()
    log_message("Outputting sorted and merged bed file")
    if not os.path.isdir("genes"):
        os.mkdir("genes")
    bed_df.to_csv(f"genes/{family}.bed", sep="\t", header=False, index=False)
    os.remove("genes/temp.bed")


if __name__ == "__main__":
    family = snakemake.wildcards.family
    hpo = snakemake.input.hpo
    ensembl = snakemake.input.ensembl
    refseq = snakemake.input.refseq
    hgnc = snakemake.input.hgnc
    main(family, hpo, ensembl, refseq, hgnc)
