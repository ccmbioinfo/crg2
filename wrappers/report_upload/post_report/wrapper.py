import requests
from json import loads
import pandas as pd
import os
import sys
import time
from PTQuery import *
from report_reshape_upload import *


logfile = snakemake.log[0]
logging.basicConfig(
    filename=logfile,
    filemode="w",
    level=logging.DEBUG,
    format="%(asctime)s:%(message)s",
    datefmt="%Y-%m-%d %H:%M",
)

report = str(snakemake.input.report)
logging.info(f"Report name: {report}")
ids = pd.read_csv(snakemake.input.ids)
logging.info(f"Participant IDs: {ids}")

credentials = pd.read_csv(snakemake.params.credentials)

data = {
    "client_id": credentials["client_id"][0],
    "username": credentials["username"][0],
    "password": credentials["password"][0],
    "grant_type": credentials["grant_type"][0],
    "scope": credentials["scope"][0],
}

# staging
base_url = "https://phenotips.genomics4rd.ca"
auth0_url = "https://genomics4rd-phenotips.us.auth0.com/oauth/token"

bearer_token = get_bearer_token(data, url=auth0_url)
BASE_REQUEST_ARGS = {"headers": {"authorization": "Bearer {}".format(bearer_token)}}

# establish a connection with Phenotips
logging.info(f"Connecting to Phenotips...")
query = PTQuery(
    base_url=base_url, base_request_args=BASE_REQUEST_ARGS, bearer_token=bearer_token
)

## Adpated from Delvin's code to Post report

# Directory to store the demultiplexed reports
folder_to_store_csvs = "report_upload/demultiplexed_reports/"

# Create the directory if it doesn't exist
logging.info(f"Creating report directory {folder_to_store_csvs}")
if not os.path.exists(folder_to_store_csvs):
    os.makedirs(folder_to_store_csvs, exist_ok=True)

# Dictionary to store information
report_dict = {"report_name": None, "family": None, "participants": []}
report_dict["report_name"] = report

logging.info(f"Generating normalized report")
participants, df = preprocess_report(
    report
)  # a list of participants and a normalized report

logging.info(f"Reshaping reports by participant")
for (
    family_participant_identifier,
    family,
    date_report_generated,
    sample_df,
) in reshape_reports(participants, df, report):
    ptp_dict = {}

    report_dict["family"] = family

    try:
        pt_id = ids.loc[ids["eids"] == family_participant_identifier, "pids"].values[0]
        ptp_dict["eid"] = family_participant_identifier
        ptp_dict["iid"] = pt_id
        logging.info(
            f"Phenotips identifier for {family_participant_identifier}: {pt_id}"
        )
    except IndexError:
        logging.error(
            f"No Phenotips identifier found for {family_participant_identifier}"
        )
        ptp_dict["variants_found"] = None
        ptp_dict["variants_removed"] = None
        ptp_dict["missing_cols"] = None
        ptp_dict["extra_cols"] = None
        ptp_dict["post_status_code"] = None

    if pt_id:

        # query variant-source-files/metadata endpoint to retrieve IDs of patients who have had reports uploaded already
        metadata = query.get_job_metadata_for_patient(pt_id)
        try:
            report = metadata[0].get("fileName")
            logging.info(
                f"Variants found for {family_participant_identifier} in report: {report}"
            )
            ptp_dict["variants_found"] = 1
            ptp_dict["variants_removed"] = 1
            logging.info(f"Removing variants for {family_participant_identifier}")
            query.delete_report(pt_id, report)
            # add a four minute wait so subsequent report upload does not fail with 409 error
            time.sleep(240)
        except IndexError:
            ptp_dict["variants_found"] = 0
            ptp_dict["variants_removed"] = 0
            logging.info(f"No variants found for {pt_id}")

        # format and POST report
        logging.info(f"Formatting report")
        extra_cols, missing_cols, formatted_ptp_report = format_for_phenotips(sample_df)

        ptp_dict["missing_cols"] = list(missing_cols)
        ptp_dict["extra_cols"] = list(extra_cols)

        ptp_fn = os.path.join(
            folder_to_store_csvs,
            f"{family_participant_identifier}_{date_report_generated}-formatted.csv",
        )

        formatted_ptp_report.to_csv(ptp_fn, index=False)

        logging.info(f"POST variants for {family_participant_identifier} to Phenotips")
        post_status_code = query.clean_and_post_report(pt_id, ptp_fn)
        ptp_dict["post_status_code"] = post_status_code
        report_dict["participants"].append(ptp_dict)

        if post_status_code != 200:
            message = f"Report POST failed for {family_participant_identifier} with code {post_status_code}"
            logging.info(f"Summary: {report_dict}")
            logging.error(message)
            sys.exit(message)
        else:
            message = f"Report POST successful for {family_participant_identifier} with code {post_status_code}"
            logging.info(message)
logging.info(f"Summary: {report_dict}")
