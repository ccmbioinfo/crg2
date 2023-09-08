from json import loads
import logging
import pandas as pd
import requests
import sys


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


def get_PT_patient_ID(C4R: str, base_url: str, token: str) -> str:
    """Get the phenotips internal ID
    API endpoint: https://docs.phenotips.com/reference/getpatientbyeid-1
    Input: C4R ID (family_participant), token"""
    logging.info(f"Querying Phenotips endpoint for participant {C4R}")
    pid = requests.get(
        f"{base_url}/rest/patients/eid/{C4R}",
        headers={"authorization": "Bearer {}".format(token)},
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


def get_PT_family_ID(pid: str, base_url: str, token: str) -> str:
    """Get the phenotips internal family ID
    API endpoint: https://docs.phenotips.com/reference/getpatient-1
    Input: project/family name, token"""
    logging.info(f"Querying Phenotips endpoint for participant {pid}")
    family = requests.get(
        f"{base_url}/rest/patients/{pid}/family",
        headers={"authorization": "Bearer {}".format(token)},
    )

    if family.status_code == 200:
        family_id = family.json()["id"]
        logging.info(f"Query successful; participant {pid} has family id {family_id}")
    else:
        # if participant is not associated with a family ID, exit with error
        logging.error(
            f"Query unsuccessful; participant {pid} is not associated with a family ID in Phenotips"
        )
        sys.exit()
    return family_id


def get_HPO(proband_id: str, base_url: str, token: str) -> pd.DataFrame:
    """Query G4RD phenotips to get suggested genes derived from HPO terms for the proband"""
    hpo = requests.get(
        f"{base_url}/rest/patients/{proband_id}/suggested-gene-panels",
        headers={"authorization": "Bearer {}".format(token)},
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


def get_pedigree_info(family_id: str, base_url: str, token: str) -> dict:
    """
    Get pedigree for a family given a Phenotips family ID
    API endpoint: https://docs.phenotips.com/reference/getfamilypedigree-1
    JSON response contains pedigree and family objects
    The family object contains only those individuals in the pedigree that have been assigned Phenotips IDs
    The pedigree object contains members, an array of pedigree nodes where each node represents a family member,
    and the relationship object, an array describing the relationships between pedigree nodes (parent:child relationships)
    """
    ped = requests.get(
        f"{base_url}/get/PhenoTips/FamilyPedigreeInterface?action=familyinfo&document_id={family_id}",
        headers={"authorization": "Bearer {}".format(token)},
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
        node_id = [
            member["id"]
            for member in ped["pedigree"]["members"]
            if member["properties"].get("id", None) == pid
        ][
            0
        ]  # pedigree node id
        node_to_C4R[node_id] = C4R
        # refer to members dict to retrieve affected status, if it exists
        affected = get_affected(pid, ped["pedigree"]["members"])
        # refer to members dict to retrieve sex
        sex = [
            member["properties"].get("sex", None)
            for member in ped["pedigree"]["members"]
            if member["properties"].get("id", None) == pid
        ][0]
        if not sex:
            sex = get_sex(pid, base_url, token)
        # refer to relationships dict to retrieve parent node IDs
        parents = get_parents(node_id, ped["pedigree"]["relationships"])
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

    return members


def get_sex(pid: int, base_url: str, token: str) -> str:
    """Get sex of an individual given Phenotips ID of that individual"""
    logging.info(f"Querying Phenotips endpoint for participant {pid}")
    pid_response = requests.get(
        f"{base_url}/rest/patients/{pid}",
        headers={"authorization": "Bearer {}".format(token)},
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


def get_affected(pid: int, members: dict) -> str:
    """Get affected status of individual"""
    affected = None
    for member in members:
        if member["properties"].get("id", None) == pid:
            try:
                affected = member["properties"]["clinicalStatus"]
                break
            except KeyError:
                try:
                    affected = member["pedigreeProperties"]["carrierStatus"]
                except KeyError:
                    affected = None

    return affected


def get_parents(node_id: int, relationships: dict) -> list:
    """
    Given node ID, retrieve IDs of parent nodes if they exist
    """
    for relationship in relationships:
        parents = None
        if node_id in [id["id"] for id in relationship["children"]]:
            parents = [id for id in relationship["members"]]
            break

    return parents


def write_pedigree(members: dict, C4R_family: str) -> None:
    """
    Write a pedigree text file given dictionary derived from Phenotips pedigree JSON
    """
    with open(f"{C4R_family}_pedigree.ped", "w") as f:
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
