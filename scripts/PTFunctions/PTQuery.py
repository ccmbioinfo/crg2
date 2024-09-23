from json import dumps, loads
import logging
from os import path
import pandas as pd
import requests
import sys
from typing import Optional


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
    # "seen_in_c4r_samples",
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


def get_bearer_token(
    data: dict, url="https://genomics4rd-phenotips.us.auth0.com/oauth/token"
) -> str:
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
                # auth=self.request_auth,
            )
        else:
            res = requests.get(
                f"{self.base_url}/rest/patients/{patient_id}",
                **self.base_request_args,
                # auth=self.request_auth,
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
                        "updateAnnotations": "true",
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
            return []

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
            # auth=self.request_auth,
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

    def get_PT_id(self, family: str, eid_df: pd.DataFrame) -> pd.DataFrame:
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
                f"{self.base_url}/rest/patients/eid/{eid}",
                auth=self.request_auth,
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

    def get_PT_patient_ID(self, C4R: str) -> str:
        """Get the phenotips internal ID
        API endpoint: https://docs.phenotips.com/reference/getpatientbyeid-1
        Input: C4R ID (family_participant), token"""
        logging.info(f"Querying Phenotips endpoint for participant {C4R}")
        pid = requests.get(
            f"{self.base_url}/rest/patients/eid/{C4R}",
            auth=self.request_auth,
        )
        if pid.status_code == 200:
            pid = pid.json()["id"]
            logging.info(f"Query successful; participant {C4R} has Phenotips id {pid}")
        else:
            # if any participant is not present in Phenotips, exit with error
            logging.error(
                f"Query unsuccessful; participant {C4R} is not present in Phenotips"
            )

        return pid

    def get_PT_family_ID(self, pid: str) -> str:
        """Get the phenotips internal family ID
        API endpoint: https://docs.phenotips.com/reference/getpatient-1
        Input: project/family name, token"""
        logging.info(f"Querying Phenotips endpoint for participant {pid}")
        family = requests.get(
            f"{self.base_url}/rest/patients/{pid}/family",
            auth=self.request_auth,
        )

        if family.status_code == 200:
            family_id = family.json()["id"]
            logging.info(
                f"Query successful; participant {pid} has family id {family_id}"
            )
        else:
            # if participant is not associated with a family ID, exit with error
            logging.error(
                f"Query unsuccessful; participant {pid} is not associated with a family ID in Phenotips"
            )
            sys.exit()
        return family_id

    def get_HPO(self, proband_id: str) -> pd.DataFrame:
        """Query G4RD phenotips to get suggested genes derived from HPO terms for the proband"""
        hpo = requests.get(
            f"{self.base_url}/rest/patients/{proband_id}/suggested-gene-panels",
            auth=self.request_auth,
        )

        hpo = hpo.json()

        hpo_df = pd.DataFrame(
            columns=[
                "Gene Symbol",
                "Gene ID",
                "Number of occurrences",
                "Features",
                "HPO IDs",
            ]
        )

        for row in hpo["rows"]:
            terms = row["terms"]
            name_translated = []
            name = []
            hpo_id = []
            for term in terms:
                name_translated.append(term["name_translated"])
                name.append(term["name"])
                hpo_id.append(term["id"])
            gene_symbol = row["gene_symbol"]
            count = row["count"]
            gene_id = row["gene_id"]

            hpo_df = pd.concat(
                [
                    hpo_df,
                    pd.DataFrame(
                        [
                            {
                                "Gene Symbol": gene_symbol,
                                "Gene ID": gene_id,
                                "Number of occurrences": int(count),
                                "Features": ", ".join(name_translated),
                                "HPO IDs": ", ".join(hpo_id),
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )
        if len(hpo_df) == 0:
            logging.info(f"Warning: participant {proband_id} has no HPO terms")
        return hpo_df

    def get_pedigree_info(self, family_id: str) -> [dict, str]:
        """
        Get pedigree for a family given a Phenotips family ID
        API endpoint: https://docs.phenotips.com/reference/getfamilypedigree-1
        JSON response contains pedigree and family objects
        The family object contains only those individuals in the pedigree that have been assigned Phenotips IDs
        The pedigree object contains members, an array of pedigree nodes where each node represents a family member,
        and the relationship object, an array describing the relationships between pedigree nodes (parent:child relationships)
        """
        ped = requests.get(
            f"{self.base_url}/get/PhenoTips/FamilyPedigreeInterface?action=familyinfo&document_id={family_id}",
            auth=self.request_auth,
        )
        ped = ped.json()
        members = (
            {}
        )  # create dictionary where keys are C4R IDs, and values are Phenotips ID, node ID, affected status, sex, and parental node IDs
        node_to_C4R = {}  # map pedigree node to C4R ID
        # iterate through familyMembers dict, rather than members dict, to retrieve only IDs of those with Phenotips ID (others are placeholders in pedigree)
        for member in ped["family"]["familyMembers"]:
            C4R = member["identifier"]
            pid = member["id"]  # Phenotips ID
            for member in ped["pedigree"]["members"]:
                properties = member.get("properties", None)
                if properties:
                    if properties.get("id", None) == pid:
                        node_id = member["id"] # pedigree node id
                        sex = member["properties"].get("sex", None)
                        break
            if not sex:
                sex = self.get_sex(pid)

            node_to_C4R[node_id] = C4R
            # refer to members dict to retrieve affected status, if it exists
            affected = self.get_affected(pid, ped["pedigree"]["members"])
            # refer to relationships dict to retrieve parent node IDs
            parents = self.get_parents(node_id, ped["pedigree"]["relationships"])
            members[C4R] = {
                "pid": pid,
                "node_id": node_id,
                "affected": affected,
                "sex": sex,
                "parents": parents,
            }
        # get C4R IDs for parents
        for member in members:
            parents = members[member]["parents"]
            if parents:
                parents_C4R = []
                for parent in parents:
                    try:
                        parents_C4R.append(node_to_C4R[parent])
                    except KeyError:
                        # parent is in pedigree but not assigned G4RD ID
                        pass
                members[member]["parents_C4R"] = parents_C4R
            else:
                members[member]["parents_C4R"] = None

        # get proband id
        proband_id = ped['pedigree']['proband']
        proband_id = node_to_C4R[proband_id]
        print(f"Proband in pedigree is {proband_id}")
        
        return members, proband_id

    def get_sex(self, pid: int) -> str:
        """Get sex of an individual given Phenotips ID of that individual"""
        logging.info(f"Querying Phenotips endpoint for participant {pid}")
        pid_response = requests.get(
            f"{self.base_url}/rest/patients/{pid}",
            auth=self.request_auth,
        )
        if pid_response.status_code == 200:
            sex = pid_response.json()["sex"]
            logging.info(f"Query successful; participant {pid} has sex {sex}")
        else:
            # if any participant is not present in Phenotips, exit with error
            logging.error(
                f"Query unsuccessful; participant {pid} is not present in Phenotips"
            )

        return sex

    def get_affected(self, pid: int, members: dict) -> str:
        """Get affected status of individual"""
        affected = None
        for member in members:
            properties = member.get("properties", None)
            if properties:
                if properties.get("id", None) == pid:
                    pedigree_properties = member.get("pedigreeProperties", None)
                    if pedigree_properties:
                        affected = member["pedigreeProperties"].get("carrierStatus", None)
                        if not affected:
                            affected = member["pedigreeProperties"].get("carrierStatus", None)
                        break

        return affected

    def get_parents(self, node_id: int, relationships: dict) -> list:
        """
        Given node ID, retrieve IDs of parent nodes if they exist
        """
        for relationship in relationships:
            parents = None
            if node_id in [id["id"] for id in relationship["children"]]:
                parents = [id for id in relationship["members"]]
                break

        return parents

    def write_pedigree(self, members: dict, C4R_family: str) -> None:
        """
        Write a pedigree text file given dictionary derived from Phenotips pedigree JSON
        """
        with open(f"/home/ccmmarvin/gene_data/pedigrees/{C4R_family}_pedigree.ped", "w") as f:
            for member in members:
                family_id = member.split("_")[0]
                sample_id = member
                member = members[member]
                if member["sex"] == "M":
                    sex = 1
                elif member["sex"] == "F":
                    sex = 2
                else:
                    sex = "other"
                affected = member["affected"]
                phenotype = 2 if affected == "affected" else 0
                paternal_id = 0
                maternal_id = 0
                if member["parents_C4R"]:
                    for p in member["parents_C4R"]:
                        parent_sex = members[p]["sex"]
                        if parent_sex == "M":
                            paternal_id = p
                        else:
                            maternal_id = p
                f.write(
                    f"{family_id} {sample_id} {paternal_id} {maternal_id} {sex} {phenotype}\n"
                )
