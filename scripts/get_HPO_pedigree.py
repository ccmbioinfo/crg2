import argparse
import logging
import pandas as pd
from PTFunctions import PTQueries
import sys


def main(C4R_ID: str, base_url: str, bearer_token: str, logfile: str) -> None:
    """
    Given a C4R ID (family_participant), retrieve family pedigree and proband HPO terms
    """

    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.DEBUG,
        format="%(asctime)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M",
    )


    # get pedigree
    logging.info(f"Retrieving pedigree for {C4R_ID}")
    C4R_family = C4R_ID.split("_")[0]
    pid = PTQueries.get_PT_patient_ID(C4R_ID, base_url, bearer_token)
    family = PTQueries.get_PT_family_ID(pid, base_url, bearer_token)
    pedigree = PTQueries.get_pedigree_info(family, base_url, bearer_token)
    PTQueries.write_pedigree(pedigree, C4R_family)

    # get HPO terms
    logging.info(f"Retrieving HPO terms for {C4R_ID}")
    HPO_df = PTQueries.get_HPO(pid, base_url, bearer_token)
    HPO_df.to_csv(
        f"{C4R_family}_HPO.txt",
        sep="\t",
        index=False,
    )
    # HPO_df.to_csv(
    #     f"/home/ccmmarvin/gene_data/HPO/{C4R_family}_HPO.txt", sep="\t", index=False
    # )


description = """Accepts TCAG batch tsv (e.g. /hpf/largeprojects/ccmbio/ccmmarvin_shared/tcag_downloads/6VBFF92-samples_for_crg2-new.tsv)
or a single C4R ID (family_participant).
Queries G4RD Phenotips to retrieve family pedigree and proband HPO terms for each family in the tsv, or for the supplied C4R ID. 
"""
parser = argparse.ArgumentParser(description=description)
parser.add_argument(
    "-f",
    "--file",
    type=str,
    required=False,
    help="TCAG batch tsv file",
)
parser.add_argument(
    "-i",
    "--id",
    type=str,
    required=False,
    help="C4R ID",
)
parser.add_argument(
    "-c",
    "--credentials",
    type=str,
    required=True,
    help="Path to Phenotips API credentials CSV",
)

args = parser.parse_args()
families = args.file
participant = args.id
credentials = args.credentials

credentials = pd.read_csv(credentials)
data = {
    "client_id": credentials["client_id"][0],
    "username": credentials["username"][0],
    "password": credentials["password"][0],
    "grant_type": credentials["grant_type"][0],
    "scope": credentials["scope"][0],
}

logging.info("Connecting to Phenotips...")
base_url = "https://phenotips.genomics4rd.ca"
bearer_token = PTQueries.get_bearer_token(data)

if families:
    # load genomes
    C4R_ids = pd.read_csv(
        families, names=["Family", "Sample", "R1", "R2", "Bam"], sep="\t"
    )
    C4R_ids = C4R_ids.drop_duplicates(subset="Family", keep="first")
    C4R_ids["C4R_id"] = (
        C4R_ids["Family"].astype(str) + "_" + C4R_ids["Sample"].astype(str)
    )
    for C4R_id in C4R_ids["C4R_id"].values.tolist():
        main(
            C4R_id,
            base_url,
            bearer_token,
            f"{C4R_id}_pedigree_HPO.log",
        )
elif participant:
    main(
        participant,
        base_url,
        bearer_token,
        f"{participant}_pedigree_HPO.log",
    )
else:
    print("Either a tsv or C4R ID must be supplied.")
    sys.exit()
