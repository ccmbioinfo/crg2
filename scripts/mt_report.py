import pandas as pd
import logging
import numpy as np
import io
import gzip
from datetime import datetime 


def log_message(*message):
    """write message to logfile and stdout"""
    if message:
        for i in message:
            logging.info(i)
            print(i)


def get_INFO_annot(df):
    """Return a dataframe with ClinVar and GnomAD annotations"""

    # dictionary to store the annotations we want to fetch from the INFO column
    data = {
        "clinvar_pathogenic": [],
        "gnomAD_AC_hom": [],
        "gnomAD_AC_het": [],
        "gnomAD_AF_hom": [],
        "gnomAD_AF_het": [],
        "gnomAD_max_hl": [],
    }

    for index, row in df.iterrows():
        INFO = row["INFO"].split(";")  # Access the INFO column and split by ";"
        INFO_names = [x.split("=")[0] for x in INFO]  # store the name of the info field
        INFO_values = [
            i.split("=")[1] for i in INFO
        ]  # store the value of the info field

        # If the key of data dict is present in INFO_names, then append corresponding INFO_value to the key's value,
        # else append value NA
        for key, value in data.items():
            if key in INFO_names:
                value.append(INFO_values[INFO_names.index(key)])
            else:
                value.append(".")

    # Create a df from the data dictionary
    annotations_df = pd.DataFrame(data)
    log_message(
        "Successfully fetched relevant ClinVar and GnomAD annotations from the INFO field of the report."
    )

    return annotations_df

def change_annot_9155(df):
    """Change status_mitomap column for m.9155A>G from . to Confirmed"""
    variant_entry = df.loc[df.HGVS == "m.9155A>G"]
    if variant_entry.empty:
        log_message("Variant m.9155A>G not present in the report.")
        return df
    else:
        df.loc[df.HGVS == "m.9155A>G", "status_mitomap"] = "Confirmed"
        df.loc[df.HGVS == "m.9155A>G", "disease_mitomap"] = "MIDD, renal insufficiency"
        log_message(
            "Found variant m.9155A>G in the report and updated status_mitomap to confirmed."
        )
        return df

def concat_df(df1, df2):
    """Concatenate two dataframes along axis 1 (column)"""

    concatenated_df = pd.concat([df1, df2], axis=1)
    log_message("Successfully joined the two dataframes")
    return concatenated_df

def remove_cols(df):
    """Remove unwanted columns from the report dataframe"""

    # List of column to be removed from the file
    remove_cols = [
        "TIER",
        "REF DEPTH",
        "TOTAL LOCUS DEPTH",
        "VARIANT QUALITY",
        "LOCUS MITOMAP",
        "QUAL",
    #    "FILTER",
        "MQM_INFO",
        "MQMR_INFO",
        "QA_INFO",
        "QR_INFO",
        "SAF_INFO",
        "SAR_INFO",
        "SRF_INFO",
        "SRR_INFO",
        "SBR_INFO",
        "SBA_INFO",
        "POS_FILTER",
        "SBR_FILTER",
        "SBA_FILTER",
        "MQMR_FILTER",
        "AQR_FILTER",
        "GT_FORMAT",
        "QR_FORMAT",
        "AQR_FORMAT",
        "QA_FORMAT",
        "AQA_FORMAT",
        "INFO",
        "FORMAT",
    ]

    # remove columns and store the remaining cols in new_df
    new_df = df.drop(remove_cols, axis=1)
    log_message("Removed unwanted columns from the dataframe")
    return new_df


def split_cols_by_sample(grouped_df):
    samples_info = {}
    for variant, row in grouped_df.iterrows():
        sample = str.split(row.SAMPLE, ";")
        ALT_DEPTH = str.split(row["ALT DEPTH"], ";")
        SAMPLE_DEPTH = str.split(row["TOTAL SAMPLE DEPTH"], ";")
        VARIANT_HETEROPLASMY = str.split(row["VARIANT HETEROPLASMY"], ";")

        x = {}
        for i in range(0, len(sample)):
            x[f"{sample[i]}.VARIANT HETEROPLASMY"] = VARIANT_HETEROPLASMY[i]
            x[f"{sample[i]}.ALT DEPTH"] = ALT_DEPTH[i]
            x[f"{sample[i]}.TOTAL SAMPLE DEPTH"] = SAMPLE_DEPTH[i]
        samples_info[variant] = x
    df = pd.DataFrame(samples_info).transpose().sort_index(axis=1).fillna("-")
    df["HGVS"] = df.index
    df["SAMPLE"] = grouped_df["SAMPLE"]

    return df


def sort_by_sample(df):
    subset_df = df[
        [
            "HGVS",
            "SAMPLE",
            "ALT DEPTH",
            "TOTAL SAMPLE DEPTH",
            "VARIANT HETEROPLASMY",
        ]
    ]
    grouped_df = subset_df.groupby("HGVS").agg(
        {
            "SAMPLE": lambda x: ";".join(str(i) for i in x),
            "ALT DEPTH": lambda x: ";".join(str(i) for i in x),
            "TOTAL SAMPLE DEPTH": lambda x: ";".join(str(i) for i in x),
            "VARIANT HETEROPLASMY": lambda x: ";".join(str(i) for i in x),
        }
    )

    df = df.drop(
        [
            "SAMPLE",
            "ALT DEPTH",
            "TOTAL SAMPLE DEPTH",
            "VARIANT HETEROPLASMY",
        ],
        1,
    )

    df2 = split_cols_by_sample(grouped_df)

    final = pd.merge(df, df2, on="HGVS", how="outer")
    log_message("Report sorted by samples")

    return final.drop_duplicates(ignore_index=True)

def get_vcf_info(vcf,report,samples):
    for i in samples:
        sample_depths=[]
        vafs=[]
        alt_depths=[]
        for row in report.iterrows():
            pos=row[1]["POS"]
            ref=row[1]["REF"]
            alt=row[1]["ALT"]
            #If pos, ref and alt match with the respective columns in the VCF then get the sample depth for that sample
            depth=list(vcf[(vcf["POS"]==pos) & (vcf["REF"]==ref) & (vcf["ALT"]==alt)][i])[0].split(":")[1]
            vaf=list(vcf[(vcf["POS"]==pos) & (vcf["REF"]==ref) & (vcf["ALT"]==alt)][i])[0].split(":")[9]
            AD=list(vcf[(vcf["POS"]==pos) & (vcf["REF"]==ref) & (vcf["ALT"]==alt)][i])[0].split(":")[5]
            sample_depths.append(depth)
            vafs.append(vaf)
            alt_depths.append(AD)
        report[f"{i}.TOTAL SAMPLE DEPTH"]=sample_depths
        report[f"{i}.VARIANT HETEROPLASMY"]=vafs
        report[f"{i}.ALT DEPTH"]=alt_depths

    return report

def check_sort(vcf,df):
    sample = df.SAMPLE.unique()
    if len(sample) == 1:
        log_message("Only one sample present in report")
        return df
    else:
        log_message("Multiple samples present in report")
        updated_df = sort_by_sample(df)
        updated_df=get_vcf_info(vcf,updated_df,sample)
        return updated_df

def remove_blacklist_pos(report):
    report=report[report["FILTER"]=="PASS"]
    return report

def reorder_cols(df):
    """Reorder columns in the report dataframe"""

    colnames = df.columns

    variant_heteroplasmy = [x for x in colnames if x.endswith("VARIANT HETEROPLASMY")]
    alt_depth = [x for x in colnames if x.endswith("ALT DEPTH")]
    total_sample_depth = [x for x in colnames if x.endswith("TOTAL SAMPLE DEPTH")]

    df=remove_blacklist_pos(df)

    col_list = [
        "CHR",
        "POS",
        "REF",
        "ALT",
        "HGVS",
        "GENE/LOCUS",
        "GENE/LOCUS DESCRIPTION",
        "COHORT COUNT",
        "SAMPLE",
    ]

    col_list2 = [
        "DISEASE MITOMAP",
        "STATUS MITOMAP",
        "DISEASE AMINO ACID CHANGE MITOMAP",
        "ALLELE FREQUENCY MITOMAP",
        "GENBANK FREQUENCY MITOMAP",
        "clinvar_pathogenic",
        "gnomAD_AC_hom",
        "gnomAD_AC_het",
        "gnomAD_AF_hom",
        "gnomAD_AF_het",
        "gnomAD_max_hl",
        "HOMOPLASMY MITOMAP",
        "HETEROPLASMY MITOMAP",
        "NUMBER OF REFERENCES MITOMAP",
        "VARIANT AMINO ACID CHANGE MITOMAP",
        "CODON POSITION MITOMAP",
        "CODON NUMBER MITOMAP",
        "NUM DISEASE REFERENCES MITOMAP",
        "RNA MITOMAP",
        "COMMERCIAL PANELS",
        "PHYLOTREE HAPLOTYPE",
        "MITOTIP SCORE",
        "MITOTIP INTERPRETATION",
        "MITOTIP PERCENTILE",
        "ANTICODON",
        "MGRB FREQUENCY",
        "MGRB FILTER",
        "MGRB AC",
        "MGRB AN",
        "PHYLOTREE MUT",
    ]
    final_col_list = (
        col_list + variant_heteroplasmy + alt_depth + total_sample_depth + col_list2
    )

    reordered_df = df[final_col_list]

    # replace '.'/'-' with '0' for some columns
    replace_col_values = [
        "ALLELE FREQUENCY MITOMAP",
        "GENBANK FREQUENCY MITOMAP",
        "gnomAD_AC_hom",
        "gnomAD_AC_het",
        "gnomAD_AF_hom",
        "gnomAD_AF_het",
        "gnomAD_max_hl",
        "MGRB FREQUENCY",
        #"MGRB FILTER",
        "MGRB AC",
        "MGRB AN",
    ]
    for col in replace_col_values:
        reordered_df[col] = reordered_df[col].replace(".", 0)

    log_message(
        "Replaced . and - with 0 for frequency columns and rearanged the columns in the dataframe"
    )
    return reordered_df

def read_vcf(vcf):
    with gzip.open(vcf, 'r') as f:
        lines = [l for l in f if not l.startswith(b'##')]
        str_lines=[]
        for l in lines:
            str_lines.append(l.decode())

    vcf_df=pd.read_csv(
        io.StringIO(''.join(str_lines)),
        sep='\t'
    )

    return vcf_df

def main(vcf, report, family):
    logfile = f"logs/report/mitochondrial/{family}.mitochondrial.report.log"
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M",
    )

    report_df = pd.read_excel(report,engine="openpyxl")
    vcf_df=read_vcf(vcf)

    CV_gnomad_annots_df = get_INFO_annot(report_df)
    merged_report = concat_df(report_df, CV_gnomad_annots_df)
    final_report = remove_cols(merged_report)
    final_report = check_sort(vcf_df,final_report)
    final_report = reorder_cols(final_report)
    final_report = change_annot_9155(final_report)

    final_report.to_csv(
        f"report/mitochondrial/{family}.mitochondrial.report.csv", index=False
    )

    current_date=datetime.now().strftime("%Y-%m-%d")

    final_report.to_csv(
        f"report/mitochondrial/{family}.mitochondrial.report.{current_date}.csv", index=False
    )

    log_message(
        "Final formatted report containing annotated list of mitochondrial variants created!"
    )


if __name__ == "__main__":
    family = snakemake.wildcards.family
    vcf= snakemake.input.vcf
    report=snakemake.input.report
    main(vcf,report,family)
