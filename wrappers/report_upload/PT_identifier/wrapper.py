import logging
import os
import pandas as pd
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "../../../scripts/PTFunctions"))
import PTQueries

logfile = snakemake.log[0]
logging.basicConfig(
    filename=logfile,
    filemode="w",
    level=logging.DEBUG,
    format="%(asctime)s:%(message)s",
    datefmt="%Y-%m-%d %H:%M",
)

eid_df = pd.read_csv(snakemake.input.samples_ids, sep="\t")
family = snakemake.params.family
credentials = pd.read_csv(snakemake.params.credentials)

data = {
    "client_id": credentials["client_id"][0],
    "username": credentials["username"][0],
    "password": credentials["password"][0],
    "grant_type": credentials["grant_type"][0],
    "scope": credentials["scope"][0],
}

logging.info("Connecting to Phenotips")
bearer_token = PTQueries.get_bearer_token(data)

df = PTQueries.get_PT_id(
    family, eid_df, "https://phenotips.genomics4rd.ca", bearer_token
)

# Write to file
df.to_csv("report_upload/PT_ids.txt", index=False)
