import logging
import os
import pandas as pd
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "../../../scripts/PTFunctions"))
from PTQuery import *

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

base_url = "https://phenotips.genomics4rd.ca"
auth0_url = "https://genomics4rd-phenotips.us.auth0.com/oauth/token"
bearer_token = get_bearer_token(data, url=auth0_url)
BASE_REQUEST_ARGS = {"headers": {"authorization": "Bearer {}".format(bearer_token)}}

# establish a connection with Phenotips
logging.info(f"Connecting to Phenotips...")
query = PTQuery(
    base_url=base_url, base_request_args=BASE_REQUEST_ARGS, bearer_token=bearer_token
)

df = query.get_PT_id(family, eid_df)


# Write to file
df.to_csv("report_upload/PT_ids.txt", index=False)
