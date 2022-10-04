import requests
from json import loads
import pandas as pd
import os
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

report=str(snakemake.input.report)
logging.info(f"Report name: {report}")
ids=pd.read_csv(snakemake.input.ids)
logging.info(f"Participant IDs: {ids}")

credentials=pd.read_csv(snakemake.params.credentials)

data={"client_id" :credentials["client_id"][0],
      "username":credentials["username"][0],
      "password":credentials["password"][0],
      "grant_type": credentials["grant_type"][0],
      "scope" : credentials["scope"][0]}

#staging
base_url = 'https://staging.phenotips.genomics4rd.ca'
auth0_url = "https://genomics4rd-phenotips-staging.us.auth0.com/oauth/token"

bearer_token = get_bearer_token(data,url = auth0_url)
BASE_REQUEST_ARGS = {"headers":{"authorization": "Bearer {}".format(bearer_token)}}

#establish a connection with Phenotips
logging.info(f"Connecting to Phenotips...")
query = PTQuery(base_url = base_url, 
                        base_request_args= BASE_REQUEST_ARGS, 
                       bearer_token = bearer_token)

## Adpated from Delvin's code to Post report

#Directory to store the demultiplexed reports
folder_to_store_csvs = 'report_upload/demultiplexed_reports/'

#Create the directory if it doesn't exist
logging.info(f"Creating report directory {folder_to_store_csvs}")
if not os.path.exists(folder_to_store_csvs):
    os.makedirs(folder_to_store_csvs, exist_ok=True)

#Dictionary to store information
report_dict = {"report_name": None, "family": None, "participants": []}
report_dict["report_name"] = report

logging.info(f"Generating normalized report")
participants, df = preprocess_report(report) # a list of participants and a normalized report

logging.info(f"Reshaping reports by participant")
for (family_participant_identifier, family,  date_report_generated, sample_df) in reshape_reports(participants, df, report):
    ptp_dict = {}

    report_dict["family"] = family                

    try:
        pt_id=ids.loc[ids["eids"]==family_participant_identifier,"pids"].values[0]
        ptp_dict["eid"] = family_participant_identifier
        ptp_dict["iid"] = pt_id
        logging.info(f"Phenotips identifier for {family_participant_identifier}: {pt_id}")
    except IndexError:
        pt_id = None
        logging.error(f"No Phenotips identifier found for {family_participant_identifier}")
        ptp_dict["variants_found"] = None
        ptp_dict["missing_cols"] = None
        ptp_dict["extra_cols"] = None
        ptp_dict["post_status_code"] = None
        

    if pt_id:

        #variants_exist = 0 
        variants_exist=query.get_variant_info(pt_id) #This line was initially commented but I'm using it
        #as a sanity check to make sure double reports don't get uploaded.
        print("No. of variants that exist:", variants_exist)
                    
        if variants_exist:
            logging.info(f"Variants found for {family_participant_identifier}")
            ptp_dict["variants_found"] = 1

        # no variants found so report will be formatted and POSTed
        else:
            logging.info(f"No variants found in Phenotips for {family_participant_identifier}, formatting variants")
            ptp_dict["variants_found"] = 0
            extra_cols, missing_cols, formatted_ptp_report = format_for_phenotips(sample_df)

            ptp_dict["missing_cols"] = list(missing_cols)
            ptp_dict["extra_cols"] = list(extra_cols)

            ptp_fn = os.path.join(
                folder_to_store_csvs,
                f"{family_participant_identifier}_{date_report_generated}-formatted.csv",
            )

            formatted_ptp_report.to_csv(ptp_fn, index=False)

            logging.info(f"POST variants for {family_participant_identifier} to Phenotips")
            # ---------------------------------------------------------------------------
            # START - UNCOMMENT THE CHUNK BELOW TO ATTEMPT A DEMULTIPLEXED REPORT UPLOAD 

            #post_status_code = query.clean_and_post_report(pt_id, ptp_fn)

            #if post_status_code != 200:
            #    logging.error(f"Report POST failed for {family_participant_identifier} with code {post_status_code}")
            #    print(f"Report POST failed for {family_participant_identifier} with code {post_status_code}")
            #else:
            #    logging.info(f"Report POST successfuly for {family_participant_identifier} with code {post_status_code}")
            #    print(f"Report POST successfuly for {family_participant_identifier} with code {post_status_code}")
                        
                        
            #ptp_dict["post_status_code"] = post_status_code
                        
            # END - UNCOMMENT THE CHUNK ABOVE TO ATTEMPT A DEMULTIPLEXED REPORT UPLOAD 
            # ------------------------- END --------------------------------------------
        report_dict["participants"].append(ptp_dict)
logging.info(f"Summary: {report_dict}")