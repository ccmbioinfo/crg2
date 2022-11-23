import requests
from json import loads
import logging
import pandas as pd
import sys

logfile = snakemake.log[0]
logging.basicConfig(
    filename=logfile,
    filemode="w",
    level=logging.DEBUG,
    format="%(asctime)s:%(message)s",
    datefmt="%Y-%m-%d %H:%M",
)


def get_bearer_token(
    data: dict, url="https://genomics4rd-phenotips.us.auth0.com/oauth/token"
) -> str:
    """Get the bearer token"""
    headers = {"content-type": "application/x-www-form-urlencoded"}
    res = requests.post(url, data=data, headers=headers)
    res2 = loads(res.text)
    token = res2.get("id_token")
    return token


def get_PT_id(family, eid_df, base_url, token):
    """
    Get the phenotips internal ID
    API endpoint: https://docs.phenotips.com/reference/getpatientbyeid-1
    Input: project/family name, samples.tsv, base_url, token
    """
    samples = list(eid_df["sample"])

    # Create empty lists to store info external id (eid), Phenotips internal id (pid), and status code
    id_dict = {"eids": [], "pids": [], "status": []}
    for s in samples:
        # Fetch PT internal identifier by searching for family_sampleID
        eid = f"{family}_{s}"
        logging.info(f"Querying Phenotips endpoint for participant {eid}")
        pid = requests.get(
            f"{base_url}/rest/patients/eid/{eid}",
            headers={"authorization": "Bearer {}".format(token)},
        )
        if pid.status_code == 200:
            pid_id = pid.json()["id"]
            logging.info(
                f"Query successful; participant {eid} has Phenotips id {pid_id}"
            )
            id_dict["pids"].append(pid_id)
        else:
            # if any participant is not present in Phenotips, exit with error
            message = (
                f"Query unsuccessful; participant {eid} is not present in Phenotips"
            )
            logging.error(message)
            sys.exit(message)
        id_dict["status"].append(pid.status_code)
        id_dict["eids"].append(eid)

    df = pd.DataFrame.from_dict(id_dict, orient="columns")

    return df


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
bearer_token = get_bearer_token(data)

df = get_PT_id(family, eid_df, "https://phenotips.genomics4rd.ca", bearer_token)

# Write to file
df.to_csv("report_upload/PT_ids.txt", index=False)
