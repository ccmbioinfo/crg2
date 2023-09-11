import os
import logging
import numpy as np

from datetime import datetime
from glob import glob
from PTQueries import *
from typing import Tuple, List, Iterator


def read_report_csv(report: str) -> pd.DataFrame:
    """
    read in a WES report and return a pandas dataframe
    """

    if report.endswith(".csv"):
        sep = ","
    elif report.endswith(".tsv"):
        sep = "\t"
    else:
        logging.error(f"Unknown fileformat and delimiter for {report}")
        raise ValueError(f"Unknown fileformat and delimiter for {report}")

    try:
        df = pd.read_csv(report, sep=sep)
    except UnicodeDecodeError:
        df = pd.read_csv(report, encoding="latin-1", sep=sep)
    except Exception as e:
        logging.error(f"Could not read in {report}")
        raise ValueError(
            f"Report '{report}' could not be read in, please double check this is a valid csv!"
        )

    return df


def preprocess_report(report: str) -> Tuple[List[str], pd.DataFrame]:
    """
    given a report path, reads in as a dataframe and performs wrangling to normalize the dataframe
    this function was originally made to normalize the report dataframes such that they could be inserted into Stager's variant store database
    """

    df = read_report_csv(report)

    # these columns are 'wide' wrt the variants
    d = {"Zygosity": [], "Burden": [], "Alt_depths": []}

    # get all columns ending with Zygosity, Burden, or Alt_depths - these are wrt to each participant
    for key in d:
        d[key].extend([col for col in df.columns if col.startswith(key)])

    # get sample names - preserved order for genotype and trio coverage
    samples = [col.replace("Zygosity.", "").strip() for col in d["Zygosity"]]

    df.columns = map(str.lower, df.columns)

    # rename to match template
    df = df.rename(columns={"omim_gene_description": "omim_phenotype"})

    # convert variation values to lowercase
    df["variation"] = df["variation"].str.lower()

    # make None's consistent (doesn't account for 0's though)
    df = df.replace({"None": None, np.nan: None})

    # replace '0's
    for bad_col in ["omim_phenotype", "orphanet"]:
        if bad_col in df.columns:
            df[bad_col] = df[bad_col].replace({"0": None})

    for col in ["conserved_in_20_mammals", "vest3_score", "revel_score", "gerp_score"]:
        if col in df.columns:
            df[col] = df[col].replace({".": None})

    df["depth"] = df["depth"].fillna(0)

    # weird column in that there are only Yes or NA's
    if "pseudoautosomal" in df.columns:
        df["pseudoautosomal"] = df["pseudoautosomal"].replace(
            {"Yes": 1, None: 0, np.nan: 0}
        )

    # remove hyperlink formatting
    for link_col in ["ucsc_link", "gnomad_link"]:
        if link_col in df.columns:
            # extract odd values b/w quotes
            try:
                df[link_col] = df[link_col].apply(lambda x: x.split('"')[1::2][0])

            # some columns don't actually have a link and therefore aren't splittable, eg. # results/14x/1473/1473.wes.2018-09-25.csv
            except IndexError as e:
                pass

    # take the last element after splitting on ';', won't affect non ';' delimited values
    df["clinvar"] = df["clinvar"].map(lambda x: str(x).split(";")[-1])

    # if coding persists. replace with appropriate value - credit to Madeline C.
    df["clinvar"] = df["clinvar"].replace(
        {
            "255": "other",
            "0": "uncertain",
            "1": "not-provided",
            "2": "benign",
            "3": "likely-benign",
            "4": "likely-pathogenic",
            "5": "pathogenic",
            "6": "drug-response",
            "7": "histocompatability",
        },
        regex=True,
    )
    # replace all forward slashes with '|' to ensure consistency
    df["clinvar"] = df["clinvar"].replace({"\\/": "|"}, regex=True)
    df["clinvar"] = df["clinvar"].replace({"None": None, np.nan: None}).str.lower()

    # these should be ints/null
    if "number_of_callers" in df.columns:
        if "Number_of_callers" in df["number_of_callers"].astype(str).values:
            df["number_of_callers"] = None

    # these should be ints/null, see. 1615.wes.regular.2020-08-14.csv
    if "frequency_in_c4r" in df.columns:
        if "Frequency_in_C4R" in df["frequency_in_c4r"].astype(str).values:
            df["frequency_in_c4r"] = None

    if "c4r_wes_counts" in df.columns:
        df["frequency_in_c4r"] = df["c4r_wes_counts"]

    # the former will eventually be deleted from the PT schema (meaning the report will fail if the col is present)
    for col in ["seen_in_c4r_samples", "c4r_wes_samples"]:
        if col in df.columns:
            df[col] = None

    # looks like report formatting changed around 2021-01. before, only spliceai_score is present and contains the impact delimited by '|'
    # afterwards, spliceai_score contains the float and spliceai_impact contains the | delimited score
    if ("spliceai_score" in df.columns) and ("spliceai_impact" not in df.columns):
        df["spliceai_impact"] = df["spliceai_score"]
        df["spliceai_score"] = None

    # check for duplicate variants wrt to pos, ref, and alt, these seem to affect a fraction of reports eg. 1666, 1743, and 1516
    dupes = df[df.duplicated(["position", "ref", "alt"], keep=False)].shape[0]

    if dupes:
        logging.info(f"Duplicate variants found for {report}")
        logging.info(f"# of variants before: {df.shape[0]}")
        df = df.drop_duplicates(["position", "ref", "alt"])
        logging.info(f"# of variants after: {df.shape[0]}")

    return samples, df


def reshape_reports(
    samples: List[str], df: pd.DataFrame, report_fn: str
) -> Iterator[Tuple[str, str, str, pd.DataFrame]]:
    """
    creates participant-wise dataframes from a cre report
    may want to move this elsewhere since it could also be called for Stager report POSTing
    """

    for i, sample in enumerate(samples):
        sample_df = df.copy(deep=True)

        # create a dictionary corresponding to the specific samples genotype columns
        genotype_cols = {
            gt: "{}.{}".format(gt, sample).lower()
            for gt in ["zygosity", "burden", "alt_depths"]
        }

        # remove existing columns in the sample df corresponding to genotype columns
        sample_df.drop(
            list(df.filter(regex="zygosity|burden|alt_depths|trio_coverage")),
            axis=1,
            inplace=True,
        )

        # iterate over each field, creating a placeholder in participant-wise dataframe, fetching the
        # appropriate participant-genotype column from the original dataframe and assigning it to a stand alone column in the participant-wise dataframe
        for gt_field, sample_specific_field in genotype_cols.items():
            # create placeholder in sample df corresponding to zygosity, burden, and alt_depths per participant
            sample_df[gt_field] = np.nan

            if genotype_cols[gt_field] in df.columns:
                # eg. sample_df[Zygosity] = df[zygosity.1389_ch0200]
                sample_df[gt_field] = df[sample_specific_field]

        # retrieve appropriate GT location
        if "gts" in df.columns:
            gts = df["gts"].str.split(",", expand=True)

            if gts.shape[1] != len(samples):
                raise IndexError(
                    "The number of extracted genotypes does not match the number of samples"
                )

            # assume the order of the genotypes is the same as the samples by indexing and replace 'gts' ;misnomer
            sample_df["gts"] = gts[i]

        if "trio_coverage" in df.columns:
            trio_coverage = df["trio_coverage"].astype(str)

            # replace dates, eg. 2003-05-06, 2004-10-03..... these are converted to dates by excel from eg. 03-05-06, 04-10-03
            # checks for values in this column that also fit eg. 11/34/54, that way we know the former was improperly converted to dates
            if any(
                trio_coverage.str.contains(
                    "[0-9]{4}-[0-9]{2}-[0-9]{2}", regex=True, na=False
                )
            ) and any(trio_coverage.str.contains("\\/")):
                trio_coverage = trio_coverage.apply(
                    lambda x: x[2:] if x.startswith("20") and "-" in x else x
                )
                trio_coverage = trio_coverage.str.replace("-", "/")

            if any(trio_coverage.str.contains("_")):
                trio_coverage = trio_coverage.str.split("_", expand=True)

                if trio_coverage.shape[1] != len(samples):
                    raise IndexError(
                        "The number of extracted coverage fields does not match the number of samples"
                    )

                # assume the order of the genotypes is the same as the samples by indexing and replace 'trio_coverage', misnomer
                sample_df["trio_coverage"] = trio_coverage[i].astype(int)

            elif any(trio_coverage.str.contains("\\/")):
                trio_coverage = trio_coverage.str.split("\\/", expand=True)

                if trio_coverage.shape[1] != len(samples):
                    raise IndexError(
                        "The number of extracted coverage fields does not match the number of samples"
                    )

                sample_df["trio_coverage"] = trio_coverage[i].astype(int)

            else:
                # singleton
                sample_df["trio_coverage"] = trio_coverage

        # filter out homozygous reference variants or insufficent coverage variants from frequency
        sample_df = sample_df[
            (sample_df["zygosity"] != "-")
            & (sample_df["zygosity"] != "Insufficient coverage")
            & (sample_df["zygosity"] != "Insufficient_coverage")
        ]

        # replace x/y/mt invalid zygosity values
        sample_df.loc[
            ~sample_df["zygosity"].isin(["Het", "Hom", None]), "zygosity"
        ] = None

        # remove MT variants
        sample_df = sample_df[~sample_df["position"].astype(str).str.startswith("MT")]

        # add identifying information based on the report
        family = os.path.basename(report_fn).split(".")[0]

        # 'gts' and 'trio_coverage' should be 'gt' and 'coverage', respectively but are retained to be in line with the example PT singleton report
        cols_to_move = [
            "position",
            "ref",
            "alt",
            "zygosity",
            "burden",
            "alt_depths",
            "gts",
            "trio_coverage",
        ]

        sample_df = sample_df[
            cols_to_move + [col for col in sample_df.columns if col not in cols_to_move]
        ]

        date_report_generated = os.path.basename(report_fn).split(".")[-2]

        yield sample, family, date_report_generated, sample_df


def format_for_phenotips(report: pd.DataFrame) -> pd.DataFrame:
    """
    Formatting for phenotips, specifically presence of columns and adding a header
    Credits to Conor Klamann for original code.
    """

    report.rename({"position": "#position"}, axis=1, inplace=True)

    missing_cols = set(ALLOWED_FIELDS).difference(set(report.columns.values))

    for col in missing_cols:
        report[col] = None

    extra_cols = set(report.columns.values).difference(set(ALLOWED_FIELDS))

    report.drop(columns=extra_cols, inplace=True)

    return extra_cols, missing_cols, report
