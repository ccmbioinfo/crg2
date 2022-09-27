import requests
from json import loads
import pandas as pd

#def get_bearer_token(data: dict, url = "https://genomics4rd-phenotips.us.auth0.com/oauth/token") -> str: #prod
def get_bearer_token(data: dict, url = "https://genomics4rd-phenotips-staging.us.auth0.com/oauth/token") -> str:
    ''' Get the bearer token '''
    headers = {"content-type": "application/x-www-form-urlencoded"}
    res = requests.post(url, data=data, headers=headers)
    res2 = loads(res.text)
    token = res2.get("id_token") 
    return token

def get_PT_id(family,eid_df,base_url,token):
    ''' 
    Get the phenotips internal ID 
    API endpoint: https://docs.phenotips.com/reference/getpatientbyeid-1 
    Input: project/family name, samples.tsv, base_url, token
    '''
    samples=list(eid_df["sample"])

    #Create empty lists to store info external id, internal id and status code
    eids=list()
    iids=list()
    status=list()
    for s in samples:
        #Fetch PT internal identifier by searching for family_sampleID
        eid=f"{family}_{s}"
        pid=requests.get(f"{base_url}/rest/patients/eid/{eid}",
                        headers={"authorization":"Bearer {}".format(token)})
        if pid.status_code == 200:
            iids.append(pid.json()["id"])
        else:
            iids.append("")
        status.append(pid.status_code)
        eids.append(eid)
    
    df=pd.DataFrame()
    df["eids"]=eids
    df["iids"]=iids
    df["status"]=status

    return df

eid_df=pd.read_csv(snakemake.input.samples_ids, sep="\t")
family=snakemake.params.family
credentials=pd.read_csv(snakemake.params.credentials)

data={"client_id" :credentials["client_id"][0],
      "username":credentials["username"][0],
      "password":credentials["password"][0],
      "grant_type": credentials["grant_type"][0],
      "scope" : credentials["scope"][0]}

bearer_token = get_bearer_token(data)

#df=get_PT_id(family,eid_df,"https://phenotips.genomics4rd.ca",bearer_token) #prod
df=get_PT_id(family,eid_df,"https://staging.phenotips.genomics4rd.ca/",bearer_token)

#Write to file
df.to_csv("report_upload/PT_ids.txt",index=False)