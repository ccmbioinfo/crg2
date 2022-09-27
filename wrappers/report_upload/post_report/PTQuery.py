from json import dumps, loads
from os import path
from re import sub
from typing import Optional

import requests
import pandas as pd

""" headers, in order, from template singleton report """
ALLOWED_FIELDS = [
    "#position",
    "ucsc_link",
    "gnomad_link",
    "ref",
    "alt",
    "zygosity",
    "gene",
    "burden",
    "gts",
    "variation",
    "info",
    "refseq_change",
    "depth",
    "quality",
    "alt_depths",
    "trio_coverage",
    "ensembl_gene_id",
    "gene_description",
    "omim_phenotype",
    "omim_inheritance",
    "orphanet",
    "clinvar",
    "frequency_in_c4r",
    #"seen_in_c4r_samples",
    "hgmd_id",
    "hgmd_gene",
    "hgmd_tag",
    "hgmd_ref",
    "gnomad_af_popmax",
    "gnomad_af",
    "gnomad_ac",
    "gnomad_hom",
    "ensembl_transcript_id",
    "aa_position",
    "exon",
    "protein_domains",
    "rsids",
    "gnomad_oe_lof_score",
    "gnomad_oe_mis_score",
    "exac_pli_score",
    "exac_prec_score",
    "exac_pnull_score",
    "conserved_in_20_mammals",
    "spliceai_impact",
    "spliceai_score",
    "sift_score",
    "polyphen_score",
    "cadd_score",
    "vest3_score",
    "revel_score",
    "gerp_score",
    "imprinting_status",
    "imprinting_expressed_allele",
    "pseudoautosomal",
    "number_of_callers",
    "old_multiallelic",
    "uce_100bp",
    "uce_200bp",
]


class BearerAuth(requests.auth.AuthBase):
    def __init__(self, token):
        self.token = token

    def __call__(self, r):
        r.headers["authorization"] = "Bearer " + self.token
        return r


#def get_bearer_token(data: dict, url = "https://genomics4rd-phenotips.us.auth0.com/oauth/token") -> str:
def get_bearer_token(data: dict, url = "https://genomics4rd-phenotips-staging.us.auth0.com/oauth/token") -> str:
    """
    get bearer token from g4rd auth0 instance
    see: https://auth0.com/docs/authorization/flows/call-your-api-using-resource-owner-password-flow for appropriate shape of 'data'
    """

    headers = {"content-type": "application/x-www-form-urlencoded"}
    res = requests.post(url, data=data, headers=headers)

    res2 = loads(res.text)
    
    token = res2.get("id_token") 
    return token


class PTQuery:
    """
    Simple class for making requests to the PhenoTips API

    Attributes:
        base_url: The url of the PT endpoint, should end with top-level domain, no slash. E.g., phenotips.example.ca
        base_request_args: dict of kwargs to pass to the requests library. Should at minimum include {"headers": {"Authorization": <...>}}.
        username/password: Only used for staging instance as basic auth is used.
        bearer_token: Only used for production; the bearer token from auth0.

    """

    def __init__(
        self,
        base_url: str,
        base_request_args: dict,
        username: Optional[str] = None,
        password: Optional[str] = None,
        bearer_token: Optional[str] = None,
    ):
        self.base_url = base_url
        self.base_request_args = base_request_args

        if username and password and bearer_token:
            print(
                "Please provide one of 1) username and password, or 2) bearer token, not both!"
            )
            sys.exit(1)
        elif username and password:
            print("Using Basic auth")
            self.request_auth = (username, password)
        elif bearer_token:
            print("Using auth0")
            self.request_auth = BearerAuth(bearer_token)

    def get_patient_info(self, patient_id: Optional[str] = None, number=5000):
        """fetch a patient's info given their id, else fetch all patients"""
        # fetch all patient info
        if not patient_id:
            res = requests.get(
                f"{self.base_url}/rest/patients?start=0&number={number}",
                **self.base_request_args,
                #auth=self.request_auth,
            )
        else:
            res = requests.get(
                f"{self.base_url}/rest/patients/{patient_id}",
                **self.base_request_args,
                #auth=self.request_auth,
            )
        if res.ok:
            return res.json()
        else:
            return res.status_code

    def clean_and_post_report(self, patient_id: str, report_path: str) -> str:
        """clean and post report for a patient, prints response status code"""
        filename = path.basename(report_path)
        args = {
            "data": {
                "metadata": dumps(
                    {
                        "patientId": patient_id,
                        "refGenome": "GRCh37",
                        "fileName": filename,
                    }
                ),
            },
            "files": {"fileStream": (None, open(report_path, "rb"))},
            "timeout": 100,
            **self.base_request_args,
        }
        res = requests.put(
            f"{self.base_url}/rest/variant-source-files/patients/{patient_id}/files/{filename}",
            **args,
            auth=self.request_auth,
        )
        return res.status_code

    def create_patient(
        self, external_id: Optional[str] = None, body: Optional[dict] = None
    ) -> str:
        """create a patient given an eid"""

        if not external_id and body:
            """assumes the externalid is passed in"""
            if "external_id" not in body:
                return "JSON body passed in but external ID not found"
        elif external_id and not body:
            """barebones participant with just the external id"""
            body = {"external_id": external_id}

        kwargs = {
            **self.base_request_args,
            "data": dumps(body),
            "headers": {
                "Content-Type": "application/json",
                "Accept": "application/json",
                **(
                    {self.base_request_args.get("headers")}
                    if self.base_request_args.get("headers")
                    else {}
                ),
            },
        }

        res = requests.post(
            f"{self.base_url}/rest/patients/", **kwargs, auth=self.request_auth
        )
        return res.status_code

    def delete_patient(self, patient_id: str) -> str:
        res = requests.delete(
            f"{self.base_url}/rest/patients/{patient_id}",
            **self.base_request_args,
            auth=self.request_auth,
        )

        return res.status_code

    def delete_report(self, patient_id: str, filename: str) -> str:
        """remove a report from a patient record (usually b/c of error when uploading)"""
        res = requests.delete(
            f"{self.base_url}/rest/variant-source-files/patients/{patient_id}/files/{filename}",
            auth=self.request_auth,
            **self.base_request_args,
        )
        return res.status_code

    def get_job_metadata_for_patient(
        self, patient_id: str, complete=False, params={}
    ) -> str:
        """search metadata for participant's latest upload attempt"""

        PAGE_SIZE = 25

        if not params:
            params = {"patientLimit": PAGE_SIZE, "patientOffset": 0}
        if complete:
            params["procStatus"] = "COMPLETE"

        res = requests.get(
            f"{self.base_url}/rest/variant-source-files/metadata",
            params=params,
            **self.base_request_args,
            auth=self.request_auth,
        )

        returned_record_count = res.json()["meta"]["returned"]

        content = res.json().get("data", [])

        found = [
            record
            for record in content
            if record["attributes"]["patientId"] == patient_id
        ]

        if found:
            return [found["attributes"] for found in found]
        elif returned_record_count > 0:
            params["patientLimit"] += PAGE_SIZE
            params["patientOffset"] += PAGE_SIZE
            return self.get_job_metadata_for_patient(patient_id, complete, params)
        else:
            raise Exception("NOT FOUND!")

    '''def get_all_job_metadata(self, params={}) -> list[dict]:
        """paginates through all job metadata, saves and returns as a list of dictionaries. useful for sanity checking whether reports were successfully processed."""

        PAGE_SIZE = 25
        FLAG = 0
        job_metadata = []
        if not params:
            params = {"patientLimit": PAGE_SIZE, "patientOffset": 0}

        while FLAG == 0:

            res = requests.get(
                f"{self.base_url}/rest/variant-source-files/metadata",
                params=params,
                **self.base_request_args,
                #auth=self.request_auth,
            )

            total_count = res.json()["meta"]["total"]

            dat = res.json().get("data", [])

            job_metadata.extend(dat)

            if len(job_metadata) < total_count:
                params["patientLimit"] += PAGE_SIZE
                params["patientOffset"] += PAGE_SIZE
            else:
                FLAG = 1

        return job_metadata''' #was giving a def get_all_job_metadata(self, params={}) -> list[dict]: TypeError: 'type' object is not subscriptable
        # so commented this function for now

    def get_patient_external_id_by_internal_id(self, patient_id: str) -> str:
        """get the patient's external ID from the internal ID"""
        res = requests.get(
            f"{self.base_url}/rest/patients/{patient_id}",
            **self.base_request_args,
            auth=self.request_auth,
        )
        if res.ok:
            return res.json()["external_id"]
        else:
            return res.status_code

    def get_internal_id_by_external_id(self, eid: str) -> str:
        """get patient's internal ID from external ID"""
        res = requests.get(
            f"{self.base_url}/rest/patients/eid/{eid}",
            **self.base_request_args,
            auth=self.request_auth,
        )
        if res.ok:
            return res.json()["id"]
        else:
            return res.status_code

    def get_variant_count(self, max=5000) -> str:
        """return count of all variants in the store up to max"""
        params = {"limit": max}
        res = requests.get(
            f"{self.base_url}/rest/variants",
            params=params,
            **self.base_request_args,
            auth=self.request_auth,
        )
        if res.ok:
            return res.json()["meta"]["returned"]
        else:
            return res.json().get("status_code", "unknown error")

    def get_match(self, gene_name: str, ensembl_id: Optional[str] = None):
        """
        fetch a collection of matches for a given query
        note that ensemblID is not required and won't have an affect on results
        """
        kwargs = {
            **self.base_request_args,
            "data": dumps(
                {
                    "gene": {"geneName": gene_name, "ensemblID": ensembl_id},
                    "variant": {"maxFrequency": 0.05, "assemblyId": "GRCh37"},
                }
            ),
            "headers": {
                "Content-Type": "application/json",
                "Accept": "application/json",
                **self.base_request_args["headers"],
            },
        }
        res = requests.post(
            f"{self.base_url}/rest/variants/match", **kwargs, auth=self.request_auth
        )
        return res.json()

    def get_variant_info(self, patient_id: str, params={}) -> str:
        """
        fetches a collection of variants for a given internal patient ID
        """

        if not params:

            params = {
                "offset": 0,
                "limit": 10,
                "sort": "chrom::asc",
                "sort": "pos::asc",
                "filter": f"patient_ids::=::{patient_id}",
            }

        res = requests.get(
            f"{self.base_url}/rest/variants/",
            params=params,
            **self.base_request_args,
            #auth=self.request_auth,
        )
        content = res.json().get("data", [])

        returned_record_count = res.json()["meta"]["returned"]

        if returned_record_count > 0:
            found = [
                record
                for record in content
                if patient_id in record["attributes"]["patient_ids"]
            ]

            return [found["attributes"] for found in found]
        else:
            return None
