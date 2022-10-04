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
    log:
        "logs/report_upload/PT_identifer.log"
    wrapper:
        get_wrapper_path("report_upload", "PT_identifier")

#Post report
rule post_report:
    input:
        ids="report_upload/PT_ids.txt",
        report=glob.glob("report/coding/*/*.wes.regular.*.csv")
    params:
        family=config["run"]["project"],
        credentials=config["run"]["PT_credentials"]
    output:
        directory("report_upload/demultiplexed_reports")
    log:
        "logs/report_upload/post_report.log"
    wrapper:
        get_wrapper_path("report_upload", "post_report")