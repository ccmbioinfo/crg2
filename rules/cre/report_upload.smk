import glob

#Find the phenotips identifier for a sample
rule PT_identifier:
    input:
        samples_ids=config["run"]["samples"]
    params:
        family=config["run"]["project"],
        credentials=config["run"]["PT_credentials"]
    output:
        "report_upload/PT_ids.txt"
    wrapper:
        get_wrapper_path("report_upload", "PT_identifier")