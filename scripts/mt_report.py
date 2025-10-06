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
        "QUAL",
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

def keep_only_pass(report):
    report=report[report["FILTER"]=="PASS"]
    return report

def reorder_cols(df):
    """Reorder columns in the report dataframe"""

    colnames = df.columns

    variant_heteroplasmy = [x for x in colnames if x.endswith("VARIANT HETEROPLASMY")]
    alt_depth = [x for x in colnames if x.endswith("ALT DEPTH")]
    total_sample_depth = [x for x in colnames if x.endswith("TOTAL SAMPLE DEPTH")]

    df=keep_only_pass(df)

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
        "MITOMAP DISEASE AC",
        "MITOMAP DISEASE AF",
        "MITOMAP DISEASE AACHANGE",
        "MITOMAP DISEASE HOMOPLASMY",
        "MITOMAP DISEASE HETEROPLASMY",
        "MITOMAP DISEASE PUBMED IDS",
        "MITOMAP DISEASE DISEASE",
        "MITOMAP DISEASE DISEASE STATUS",
        "MITOMAP DISEASE HGFL",
        "MITOMAP CONFIRMED MUTATIONS LOCUS",
        "MITOMAP CONFIRMED MUTATIONS LOCUSTYPE",
        "MITOMAP CONFIRMED MUTATIONS ASSOCIATEDDISEASE",
        "MITOMAP CONFIRMED MUTATIONS ALLELE",
        "MITOMAP CONFIRMED MUTATIONS AMINOACIDCHANGE",
        "MITOMAP CONFIRMED MUTATIONS STATUSMITOMAPCLINGEN",
        "MITOMAP CONFIRMED MUTATIONS LASTUPDATE",
        "MITOMAP MUTATIONS CODING CONTROL LOCUS",
        "MITOMAP MUTATIONS CODING CONTROL ALLELE",
        "MITOMAP MUTATIONS CODING CONTROL DISEASE",
        "MITOMAP MUTATIONS CODING CONTROL NUCLEOTIDECHANGE",
        "MITOMAP MUTATIONS CODING CONTROL AMINOACIDCHANGE",
        "MITOMAP MUTATIONS CODING CONTROL PLASMY",
        "MITOMAP MUTATIONS CODING CONTROL STATUS",
        "MITOMAP MUTATIONS CODING CONTROL GB FREQ",
        "MITOMAP MUTATIONS CODING CONTROL GB SEQS",
        "MITOMAP MUTATIONS CODING CONTROL REFERENCES",
        "MITOMAP MUTATIONS RNA LOCUS",
        "MITOMAP MUTATIONS RNA DISEASE",
        "MITOMAP MUTATIONS RNA ALLELE",
        "MITOMAP MUTATIONS RNA RNA",
        "MITOMAP MUTATIONS RNA HOMOPLASMY",
        "MITOMAP MUTATIONS RNA HETEROPLASMY",
        "MITOMAP MUTATIONS RNA STATUS",
        "MITOMAP MUTATIONS RNA MITOTIP",
        "MITOMAP MUTATIONS RNA GB FREQ",
        "MITOMAP MUTATIONS RNA GB SEQS",
        "MITOMAP MUTATIONS RNA REFERENCES",
        "MITOMAP POLYMORPHISMS AC",
        "MITOMAP POLYMORPHISMS AF",
        "MITOMAP POLYMORPHISMS HGFL",
        "MITOMAP VARIANTS CODING LOCUS",
        "MITOMAP VARIANTS CODING NUCLEOTIDECHANGE",
        "MITOMAP VARIANTS CODING CODONNUMBER",
        "MITOMAP VARIANTS CODING CODONPOSITION",
        "MITOMAP VARIANTS CODING AMINOACIDCHANGE",
        "MITOMAP VARIANTS CODING GB FREQ",
        "MITOMAP VARIANTS CODING GB SEQS",
        "MITOMAP VARIANTS CODING CURATEDREFERENCES",
        "MITOMAP VARIANTS CONTROL LOCUS",
        "MITOMAP VARIANTS CONTROL NUCLEOTIDECHANGE",
        "MITOMAP VARIANTS CONTROL GB FREQ",
        "MITOMAP VARIANTS CONTROL GB SEQS",
        "MITOMAP VARIANTS CONTROL CURATEDREFERENCES",
        "clinvar_pathogenic",
        "gnomAD_AC_hom",
        "gnomAD_AC_het",
        "gnomAD_AF_hom",
        "gnomAD_AF_het",
        "gnomAD_max_hl",
        "PHYLOTREE HAPLOTYPE",
        "MITOTIP SCORE",
        "MITOTIP PERCENTILE",
        "MITOTIP QUARTILE",
        "MITOTIP SCORE INTERPRETATION",
        "MITOMAP STATUS",
        "COUNT",
        "PERCENTAGE",
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
        "MITOMAP POLYMORPHISMS AF",
        "MITOMAP POLYMORPHISMS AC",
        "gnomAD_AC_hom",
        "gnomAD_AC_het",
        "gnomAD_AF_hom",
        "gnomAD_AF_het",
        "gnomAD_max_hl",
        "MGRB FREQUENCY",
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

    final_report = remove_cols(report_df)
    final_report = check_sort(vcf_df,final_report)
    final_report = reorder_cols(final_report)

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
