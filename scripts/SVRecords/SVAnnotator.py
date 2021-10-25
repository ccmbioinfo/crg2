import numpy as np
import pandas as pd
import sqlite3
import os
import subprocess
from pybedtools import BedTool
from collections import defaultdict
from enum import Enum, auto

# original code written by Dennis Kao: https://github.com/ccmbioinfo/crg/blob/master/SVRecords/SVAnnotator.py

class SVTYPE(Enum):
    def _generate_next_value_(name, start, count, last_values):
        return name

    DEL = auto()
    DUP = auto()
    INS = auto()
    INV = auto()
    IDP = auto()
    TRA = auto()


class SVAnnotator:
    def __init__(self, exon_bed, hgmd_db, hpo, exac, omim, biomart):

        self.make_gene_ref_df(biomart)

        if os.path.isfile(hpo):
            print("Annotating genes with HPO terms")
            self.HPO = True
            self.annotate_hpo(hpo)
        else:
            print("No valid HPO file specified. Skipping annotation.")
            self.HPO = False

        print("Annotating genes with OMIM phenotypes and inheritance patterns")
        self.annotate_omim(omim)
        print("Annotating genes with ExAC transcript probabilities")
        self.annotate_exac(exac)

        # technical note: drop duplicates before setting index - doing the reverse order will drop all duplicated columns instead of keeping one copy
        self.gene_ref_df = (
            self.gene_ref_df.drop_duplicates(keep="first")
            .set_index("BioMart Ensembl Gene ID")
            .astype(str)
        )

    def set_column_values(self, df, annotation_dict, column_name):
        for interval, data in annotation_dict.items():
            df.loc[interval, column_name] = data

    def append_prefix_to_columns(self, df, prefix):
        # used to avoid collision of column names between different dataframes
        df.rename(
            columns={col: "%s %s" % (prefix, col) for col in df.columns}, inplace=True
        )

    def make_hgnc_dict(self, hgnc):
        hgnc_dict = {}
        with open(hgnc) as f:
            for line in f.readlines()[1:]:
                fields = line.strip("\n").split("\t")
                approved_symbol = fields[1].upper()
                prev_symbols = fields[4].split(", ")
                synonymous_symbols = fields[5].split(", ")

                for sym in prev_symbols:
                    hgnc_dict[sym.upper()] = approved_symbol

                for sym in synonymous_symbols:
                    hgnc_dict[sym.upper()] = approved_symbol

        self.hgnc_dict = hgnc_dict

    def gene2hgnc(self, col):
        def translate(cell):
            if isinstance(cell, str):
                if cell in self.hgnc_dict:
                    # print('individual sym match: %s %s' % (cell, self.hgnc_dict[cell]))
                    return self.hgnc_dict[cell]
                return cell
                # return self.hgnc_dict[cell] if cell in self.hgnc_dict else cell
            elif isinstance(cell, list):
                for sym in cell:
                    if sym in self.hgnc_dict:
                        # if len(cell) > 1: print('list match: %s %s' % (sym, self.hgnc_dict[sym]))
                        # else: print('single element match: %s %s' % (sym, self.hgnc_dict[sym]))
                        return self.hgnc_dict[sym]
                return cell[
                    0
                ]  # translation to hgnc gene has failed, just picked the first gene name then
            else:
                ValueError(
                    "DataFrame member is not of type str or list. %s" % str(cell)
                )

        return col.apply(translate)

    def calc_exons_spanned(self, sample_df, exon_bed):
        print("Calculating the number of exons affected by each structural variant ...")

        exon_counts = defaultdict(
            int
        )  # incrementing dict and setting value in df is faster than incrementing values in the df
        exon_ref = BedTool(exon_bed)
        sample_bedtool = BedTool(
            list(sample_df.reset_index()[["CHROM", "POS", "END", "SVTYPE"]].values)
        )

        for interval in sample_bedtool.intersect(exon_ref, wa=True):
            exon_counts[
                (
                    str(interval.chrom),
                    str(interval.start),
                    str(interval.stop),
                    str(interval[3]),
                )
            ] += 1

        count_df = pd.Series(exon_counts).to_frame().astype(str)
        print(count_df.head())
        count_df.index.names = ["CHROM", "POS", "END", "SVTYPE"]
        count_df.columns = ["EXONS_SPANNED"]

        return sample_df.join(count_df).fillna(value={"EXONS_SPANNED": 0})

    def calc_exon_boundaries(self, sample_df, exon_bed):
        print(
            "Calculating exon boundaries closest to structural variant breakpoints ..."
        )

        def find_min_distance(position, boundaries):
            # range is boundary minus position of breakpoint plus one to account for 1-based coordinates
            distance = {
                ((boundary - position) + 1): boundary
                for boundary in boundaries["value"].tolist()
            }
            min_distance = min(distance, key=abs)
            min_boundary = distance[min_distance]
            gene = boundaries[boundaries["value"] == min_boundary]["GENE"].values[0]
            return min_distance, min_boundary, gene

        exons = pd.read_csv(exon_bed, sep="\t", names=["CHROM", "POS", "END", "GENE"])
        # make single column containing all exon boundaries
        exons = pd.melt(exons, id_vars=["CHROM", "GENE"], value_vars=["POS", "END"])

        boundary_distances = defaultdict()
        boundary_distances["left"] = {"nearest_boundary": [], "nearest_distance": []}
        boundary_distances["right"] = {"nearest_boundary": [], "nearest_distance": []}

        for index, row in sample_df.iterrows():
            chr, pos, end = str(row["CHROM"]), int(row["POS"]), int(row["END"])
            print(chr, pos, end)
            boundaries = exons[exons["CHROM"] == chr]
            if len(boundaries) != 0:
                for breakpoint in "left", "right":
                    position = pos if breakpoint == "left" else end
                    min_distance, min_boundary, gene = find_min_distance(
                        position, boundaries
                    )
                    boundary_distances[breakpoint]["nearest_boundary"].append(
                        gene + "|" + str(min_boundary)
                    )
                    boundary_distances[breakpoint]["nearest_distance"].append(
                        min_distance
                    )
            else:
                for breakpoint in "left", "right":
                    # for non-canonical chromosomes or MT
                    boundary_distances[breakpoint]["nearest_boundary"].append(".")
                    boundary_distances[breakpoint]["nearest_distance"].append(".")

        sample_df["nearestLeftExonBoundary"], sample_df["nearestLeftExonDistance"] = (
            boundary_distances["left"]["nearest_boundary"],
            boundary_distances["left"]["nearest_distance"],
        )
        sample_df["nearestRightExonBoundary"], sample_df["nearestRightExonDistance"] = (
            boundary_distances["right"]["nearest_boundary"],
            boundary_distances["right"]["nearest_distance"],
        )

        return sample_df.set_index(["CHROM", "POS", "END", "SVTYPE"])

    def find_min_distance(self, position, boundaries):
        # range is boundary minus position of breakpoint plus one to account for 1-based coordinates
        distance = {
            ((boundary - position) + 1): boundary
            for boundary in boundaries["value"].tolist()
        }
        min_distance = min(distance, key=abs)
        min_boundary = distance[min_distance]
        gene = boundaries[boundaries["value"] == min_boundary]["GENE"].values[0]
        return min_distance, min_boundary, gene

    def annotate_hgmd(self, hgmd, sv_record):
        print(
            "Annotating genes with published cases of pathogenic structural variants from HGMD"
        )

        def get_hgmd_df():
            conn = sqlite3.connect(hgmd)

            gros_del = pd.read_sql_query(
                """
                SELECT del.DISEASE, del.TAG, del.DESCR, del.gene,
                printf('%s:%s:%s:%s', del.JOURNAL, del.AUTHOR, del.YEAR, del.PMID) AS JOURNAL_DETAILS,
                ALLGENES.HGNCID
                FROM GROSDEL as del 
                LEFT JOIN ALLGENES ON ALLGENES.GENE=del.GENE;
            """,
                conn,
            )

            gros_ins = pd.read_sql_query(
                """
                SELECT ins.DISEASE, ins.TAG, ins.DESCR, ins.gene, 
                printf('%s:%s:%s:%s', ins.JOURNAL, ins.AUTHOR, ins.YEAR, ins.PMID) AS JOURNAL_DETAILS,
                ALLGENES.HGNCID
                FROM GROSINS as ins 
                LEFT JOIN ALLGENES ON ALLGENES.GENE=ins.GENE
                WHERE ins.type='I';
            """,
                conn,
            )

            gros_dup = pd.read_sql_query(
                """
                SELECT ins.DISEASE, ins.TAG, ins.DESCR, ins.gene, 
                printf('%s:%s:%s:%s', ins.JOURNAL, ins.AUTHOR, ins.YEAR, ins.PMID) AS JOURNAL_DETAILS,
                ALLGENES.HGNCID
                FROM GROSINS as ins 
                LEFT JOIN ALLGENES ON ALLGENES.GENE=ins.GENE
                WHERE ins.type='D';
            """,
                conn,
            )

            conn.close()

            return gros_del, gros_ins, gros_dup

        def groupby_genes(df):
            df["hgncID"] = df["hgncID"].astype(str)
            # df['omimid'] = df['omimid'].astype(str)
            # df['gene'] = self.gene2hgnc(df['gene'])
            df = df.groupby(by="gene", as_index=False).agg(
                lambda x: "%s" % " & ".join(x)
            )
            df["hgncID"] = df["hgncID"].apply(lambda col: col.split(", ")[0])
            # df['omimid'] = df['omimid'].apply(lambda col: col.split(', ')[0])
            return df

        hgmd_sv_df = pd.DataFrame()
        gros_del, gros_ins, gros_dup = get_hgmd_df()

        gros_del = groupby_genes(gros_del)
        gros_ins = groupby_genes(gros_ins)
        gros_dup = groupby_genes(gros_dup)

        for df in [gros_del, gros_ins, gros_dup]:
            df["gene"] = df["gene"].apply(lambda symbol: symbol.upper())
            self.append_prefix_to_columns(df, "HGMD")

        gros_del["HGMD SVTYPE"] = SVTYPE.DEL.value
        gros_ins["HGMD SVTYPE"] = SVTYPE.INS.value
        gros_dup["HGMD SVTYPE"] = SVTYPE.DUP.value

        # hgmd_sv_df = hgmd_sv_df.rename(columns={'HGMD gene': 'Genes in HGMD'})
        hgmd_sv_df = pd.concat(
            [gros_del, gros_ins, gros_dup], ignore_index=True, sort=False
        )
        hgmd_sv_df["Genes in HGMD"] = hgmd_sv_df["HGMD gene"]
        hgmd_sv_df = hgmd_sv_df.set_index(keys=["HGMD gene", "HGMD SVTYPE"]).astype(str)

        return sv_record.join(
            hgmd_sv_df, on=["BioMart Associated Gene Name", "SVTYPE"], how="left"
        )

    def prioritized_annotation(self, gene_ref_df, annotation_df, matched_fields):
        matched_rows = []

        ann_match_columns = list(matched_fields.keys())
        ref_match_columns = list(matched_fields.values())
        gene_ref_df_matching_cols = gene_ref_df[ref_match_columns].drop_duplicates()

        annotation_df = annotation_df.drop_duplicates()

        # join on equivalent fields
        for ann_field, ref_field in matched_fields.items():
            matched = gene_ref_df_matching_cols.join(
                annotation_df.set_index(ann_field), on=ref_field, how="inner"
            )
            matched_rows.append(matched)

        # drop columns from annotation_df used for joining
        # for table in matched_rows:
        #     for field in ann_match_columns:
        #         try:
        #             table.drop(columns=[field, ], axis=0, inplace=True)
        #             #print("Dropped %s" % field)
        #         except KeyError:
        #             pass #tried to drop column which was joined on

        merged_df = (
            pd.concat(matched_rows, ignore_index=True, sort=False)
            .drop_duplicates()
            .set_index(ref_match_columns)
            .dropna(how="all")
        )

        # add the remaining fields to the reference dataframe
        return self.gene_ref_df.join(
            merged_df, on=ref_match_columns, how="left"
        ).drop_duplicates()

    def left_join(self, df1, df2, field1, field2):
        df2 = df2.drop_duplicates().set_index(field2).dropna(how="all")
        return (
            df1.set_index(field1)
            .join(df2, how="left")
            .drop_duplicates()
            .rename_axis(field1)
            .reset_index()
        )

    def annotate_hpo(self, hpo):
        matching_fields = {
            "HPO Gene ID": "BioMart Ensembl Gene ID",
        }

        hpo_df = pd.read_csv(hpo, sep="\t")
        hpo_df.columns = hpo_df.columns.str.strip()
        # hpo_df = hpo_df[['Gene ID', 'Gene symbol', 'Features']]
        hpo_df = hpo_df[["Gene ID", "Features"]]
        hpo_df["Features"] = hpo_df["Features"].apply(
            lambda features: features.replace("; ", ", ")
        )

        hpo_df = hpo_df.astype(str)
        self.append_prefix_to_columns(hpo_df, "HPO")
        hpo_df = self.left_join(
            hpo_df,
            self.gene_ref_df[
                ["BioMart Ensembl Gene ID", "BioMart Associated Gene Name"]
            ],
            "HPO Gene ID",
            "BioMart Ensembl Gene ID",
        )
        hpo_df = hpo_df.rename(columns={"BioMart Associated Gene Name": "Genes in HPO"})

        self.gene_ref_df = self.prioritized_annotation(
            self.gene_ref_df, hpo_df, matching_fields
        )
        # self.gene_ref_df.to_csv("hpo_ann.tsv", sep="\t")

    def annotate_omim(self, omim):
        omim_inheritance_codes = {
            "Autosomal dominant": "AD",
            "Autosomal recessive": "AR",
            "X-linked dominant": "XLD",
            "X-linked recessive": "XLR",
            "Y-linked dominant": "YLD",
            "Y-linked recessive": "YLR",
            "X-linked": "XL",
            "Y-linked": "YL",
        }
        omim_inheritance_codes = {
            key.lower(): value for key, value in omim_inheritance_codes.items()
        }  # for case insensitive text search

        def process_OMIM_phenotype(phenotype):
            """
            omim phenotype example:

            {Epilepsy, generalized, with febrile seizures plus, type 5, susceptibility to}, 613060 (3), Autosomal dominant;
            {Epilepsy, idiopathic generalized, 10}, 613060 (3), Autosomal dominant;
            {Epilepsy, juvenile myoclonic, susceptibility to}, 613060 (3), Autosomal dominant
            """

            inheritance = []

            if pd.isnull(phenotype):
                return phenotype
            else:
                for p in phenotype.split("; "):
                    multiple_inheritance = [
                        code
                        for description, code in omim_inheritance_codes.items()
                        if description.lower() in p.lower()
                    ]
                    if multiple_inheritance:
                        inheritance.append("&".join(multiple_inheritance))

                return ", ".join(inheritance)

        matching_fields = {"OMIM Ensembl Gene ID": "BioMart Ensembl Gene ID"}

        # OMIM adds comments to their CSV file. These comments start with '#' character and are present in the header and footer of the file.
        omim_df = pd.read_csv(omim, sep="\t", header=3, skipfooter=61, engine="python")
        omim_df.columns = omim_df.columns.str.replace("#", "")
        omim_df.columns = omim_df.columns.str.strip()

        omim_df = omim_df[["MIM Number", "Ensembl Gene ID", "Phenotypes"]]
        omim_df = omim_df[
            pd.notnull(omim_df["Phenotypes"])
        ]  # drop all nan phenotype columns
        omim_df["Inheritance"] = omim_df["Phenotypes"].apply(
            lambda col: process_OMIM_phenotype(col)
        )
        omim_df = (
            omim_df.astype(str)
            .groupby("Ensembl Gene ID", as_index=False)
            .agg(
                {
                    "Phenotypes": " & ".join,
                    "MIM Number": " & ".join,
                    "Inheritance": " & ".join,
                }
            )
        )
        self.append_prefix_to_columns(omim_df, "OMIM")

        omim_df = self.left_join(
            omim_df,
            self.gene_ref_df[
                ["BioMart Ensembl Gene ID", "BioMart Associated Gene Name"]
            ],
            "OMIM Ensembl Gene ID",
            "BioMart Ensembl Gene ID",
        )
        omim_df = omim_df.rename(
            columns={"BioMart Associated Gene Name": "Genes in OMIM"}
        )

        self.gene_ref_df = self.prioritized_annotation(
            self.gene_ref_df, omim_df, matching_fields
        )
        # self.gene_ref_df.to_csv("omim_ann.tsv", sep="\t")

    def annotate_exac(self, exac):
        matching_fields = {
            "ExAC gene": "BioMart Associated Gene Name",
        }

        exac_df = pd.read_csv(exac, sep="\t")
        exac_df.columns = exac_df.columns.str.strip()
        exac_df["transcript"] = exac_df["transcript"].apply(
            lambda transcript_id: transcript_id.split(".")[0]
        )
        exac_df = exac_df[["gene", "syn_z", "mis_z", "lof_z", "pLI"]]

        exac_df = exac_df.astype(str)
        self.append_prefix_to_columns(exac_df, "ExAC")

        self.gene_ref_df = self.prioritized_annotation(
            self.gene_ref_df, exac_df, matching_fields
        )
        # self.gene_ref_df.to_csv("exac_ann.tsv", sep="\t")

    def annotate_gnomad(self, gnomad, sv_record, reciprocal_overlap=0.5):
        print(
            "Annotating structural variants with those seen in gnomAD_SV based on a %f reciprocal overlap ..."
            % reciprocal_overlap
        )

        gnomad_cols = [
            "CHROM",
            "START",
            "END",
            "NAME",
            "SVTYPE",
            "AN",
            "AC",
            "AF",
            "N_HOMREF",
            "N_HET",
            "N_HOMALT",
            "FREQ_HOMREF",
            "FREQ_HET",
            "FREQ_HOMALT",
            "POPMAX_AF",
        ]
        gnomad_ann_cols = [
            "gnomAD_SVTYPE",
            "gnomAD_AN",
            "gnomAD_AC",
            "gnomAD_AF",
            "gnomAD_N_HOMREF",
            "gnomAD_N_HET",
            "gnomAD_N_HOMALT",
            "gnomAD_FREQ_HOMREF",
            "gnomAD_FREQ_HET",
            "gnomAD_FREQ_HOMALT",
            "gnomAD_POPMAX_AF",
        ]
        gnomad_df = pd.read_csv(gnomad, sep="\t", dtype="str").astype(str)
        gnomad_df.columns = gnomad_df.columns.str.replace("#", "")
        gnomad_df.columns = gnomad_df.columns.str.strip()
        gnomad_df = gnomad_df[gnomad_cols]
        gnomad_bed = BedTool(gnomad_df.itertuples(index=False))

        sample_sv = sv_record.make_ref_bedtool()

        ann_df = (
            sample_sv.intersect(
                gnomad_bed, wa=True, wb=True, F=reciprocal_overlap, f=reciprocal_overlap
            )
            .to_dataframe(
                names=[
                    "CHROM",
                    "POS",
                    "END",
                    "SVTYPE",
                    "gnomAD_CHROM",
                    "gnomAD_START",
                    "gnomAD_END",
                    "gnomAD_ID",
                ]
                + gnomad_ann_cols
            )
            .astype(str)
        )
        ann_df = ann_df.drop(ann_df[ann_df["SVTYPE"] != ann_df["gnomAD_SVTYPE"]].index)

        ann_df["gnomAD_SV"] = ann_df[
            ["gnomAD_CHROM", "gnomAD_START", "gnomAD_END"]
        ].apply(lambda x: "{}:{}-{}".format(x[0], x[1], x[2]), axis=1)
        ann_df = ann_df.drop(columns=["gnomAD_CHROM", "gnomAD_START", "gnomAD_END"])
        ann_df = ann_df.groupby(["CHROM", "POS", "END", "SVTYPE"]).agg(list)
        ann_df = ann_df[ann_df.columns].applymap(lambda cell: " & ".join(cell))

        return sv_record.df.join(ann_df)

    def annotate_counts(
        self, counts, sv_record, prefix="COUNT", reciprocal_overlap=0.5
    ):
        print(
            "Annotating structural variants with those seen in %s based on a %f reciprocal overlap ..."
            % (counts, reciprocal_overlap)
        )

        cols = ["COUNT_CHROM", "COUNT_START", "COUNT_END", "COUNT_SVTYPE", "COUNT"]

        count_df = pd.read_csv(counts, sep="\t", dtype="str").astype(str)
        count_bed = BedTool(count_df.itertuples(index=False))

        sample_sv = sv_record.make_ref_bedtool()

        ann_df = sample_sv.intersect(
            count_bed, wa=True, wb=True, F=reciprocal_overlap, f=reciprocal_overlap
        ).to_dataframe(
            names=[
                "CHROM",
                "POS",
                "END",
                "SVTYPE",
            ]
            + cols
        )
        ann_df[["CHROM", "POS", "END", "SVTYPE",]] = ann_df[
            [
                "CHROM",
                "POS",
                "END",
                "SVTYPE",
            ]
        ].astype(
            str
        )  # reference dataframe is typecasted as string
        ann_df = ann_df.drop(ann_df[ann_df["SVTYPE"] != ann_df["COUNT_SVTYPE"]].index)

        ann_df = ann_df.groupby(["CHROM", "POS", "END", "SVTYPE"]).agg(
            {
                "COUNT_CHROM": "first",
                "COUNT_SVTYPE": "first",
                "COUNT_START": "min",
                "COUNT_END": "max",
                "COUNT": "sum",
            }
        )
        ann_df["COUNT_SV"] = ann_df[["COUNT_CHROM", "COUNT_START", "COUNT_END"]].apply(
            lambda x: "{}:{}-{}".format(x[0], x[1], x[2]), axis=1
        )
        ann_df = ann_df.drop(
            columns=["COUNT_CHROM", "COUNT_START", "COUNT_END", "COUNT_SVTYPE"]
        )

        ann_df.columns = ann_df.columns.str.replace("COUNT", prefix)

        df = sv_record.df.join(ann_df)
        df[[prefix]] = df[[prefix]].fillna(0)

        return df

    def annotsv(self, sample_df):
        """
        Handles DGV, DDD annotations
        """
        all_sv_bed_name = "all_sv.bed"
        annotated = "./{}.annotated.tsv".format(all_sv_bed_name)
        sample_df.reset_index()[["CHROM", "POS", "END", "SVTYPE"]].to_csv(
            all_sv_bed_name, index=False, sep="\t"
        )
        subprocess.call(
            "$ANNOTSV/bin/AnnotSV -SVinputFile {} -SVinputInfo 1 -outputFile {}".format(
                all_sv_bed_name, annotated
            ),
            shell=True,
        )

        annotsv_df_original = pd.read_csv(annotated, sep="\t").astype(str)
        # DDD annotations are only given to SVs that are 'split'
        # but we want the DGV annotations that are given to the 'full' SVs
        # so, separate these and join them to get both
        sv_cols = ["SV chrom", "SV start", "SV end", "SV type"]
        DDD_cols = [
            "DDD_status",
            "DDD_mode",
            "DDD_consequence",
            "DDD_disease",
            "DDD_pmids",
        ]
        annotsv_df_split = annotsv_df_original[
            annotsv_df_original["AnnotSV type"] == "split"
        ][sv_cols + DDD_cols].set_index(sv_cols)
        annotsv_df_full = annotsv_df_original[
            annotsv_df_original["AnnotSV type"] == "full"
        ][
            sv_cols
            + [
                col
                for col in annotsv_df_original.columns.tolist()
                if "DGV" in col.upper()
            ]
        ].set_index(
            sv_cols
        )
        annotsv_df_all = (
            annotsv_df_split.merge(annotsv_df_full, how="outer", on=sv_cols)
            .reset_index()
            .drop_duplicates(subset=sv_cols)
        )
        annotsv_df_cols = [
            col for col in annotsv_df_all.columns.tolist() if col not in DDD_cols
        ]
        annotsv_df_all = annotsv_df_all.fillna(".")
        # annotsv_df_split will have multiple rows per SV (one for each gene overlapped by the SV); aggregate to de-duplicate
        annotsv_df = (
            annotsv_df_all.groupby(annotsv_df_cols)[DDD_cols]
            .agg(
                {
                    "DDD_status": ",".join,
                    "DDD_mode": ",".join,
                    "DDD_consequence": ",".join,
                    "DDD_disease": ",".join,
                    "DDD_pmids": ",".join,
                }
            )
            .reset_index()
        )

        annotsv_df = annotsv_df.rename(
            columns={
                annotsv_df.columns[0]: "CHROM",
                annotsv_df.columns[1]: "POS",
                annotsv_df.columns[2]: "END",
                annotsv_df.columns[3]: "SVTYPE",
            }
        ).set_index(keys=["CHROM", "POS", "END", "SVTYPE"])
        sample_df = sample_df.join(annotsv_df)

        # os.remove(all_sv_bed_name)
        # os.remove(annotated)

        return sample_df

    def make_gene_ref_df(self, biomart):
        df = pd.read_csv(biomart, sep="\t")
        # df = df[['Ensembl Gene ID', 'Ensembl Transcript ID', 'Associated Gene Name', 'HGNC ID(s)', 'MIM Gene Accession']].drop_duplicates()
        df = df[
            [
                "Ensembl Gene ID",
                "Associated Gene Name",
            ]
        ].drop_duplicates()
        df["Associated Gene Name"] = df["Associated Gene Name"].apply(
            lambda symbol: symbol.upper()
        )  # make all gene symbols a single case to increase match rate with other dataframes
        df = df.astype(str)
        self.append_prefix_to_columns(df, "BioMart")

        self.gene_ref_df = df

    def annotate_genes(self, sample_df, gene_col):
        def count_unique_terms(cell):
            terms = set()

            for elem in cell:
                if pd.isnull(elem):
                    continue
                elif ", " in elem:
                    terms.update(elem.split(", "))
                elif elem != "na" and elem != "nan":
                    terms.add(elem)

            return len(terms)

        # extract genes from sample_df, create a new dataframe where each row only has a single ensemble id and interval info
        gene_df = (
            sample_df.apply(lambda x: pd.Series(x[gene_col]), axis=1)
            .stack()
            .reset_index(level=4, drop=True)
        )
        gene_df = gene_df.to_frame().rename(columns={0: gene_col}).astype(str)
        # gene_df.to_csv('seperated_genes.csv')

        # annotate passed in ensemble gene id's using the generated reference dataframe
        gene_df = gene_df.join(self.gene_ref_df, on=gene_col, how="left").reset_index()
        # gene_df.to_csv('annotated_genes.csv')

        # aggregate all annotation columns within the same sv interval
        gene_df = gene_df.groupby(["CHROM", "POS", "END", "SVTYPE"]).agg(list)
        # gene_df.to_csv('grouped_index.csv')

        # add cardinality columns
        if self.HPO:
            gene_df["N_UNIQUE_HPO_TERMS"] = [
                [count_unique_terms(values["HPO Features"])]
                for index, values in gene_df.iterrows()
            ]
            gene_df["N_GENES_IN_HPO"] = [
                [count_unique_terms(values["Genes in HPO"])]
                for index, values in gene_df.iterrows()
            ]
        gene_df["N_GENES_IN_OMIM"] = [
            [count_unique_terms(values["Genes in OMIM"])]
            for index, values in gene_df.iterrows()
        ]

        # parse out and replace nan values with "na" string
        gene_df[gene_df.columns] = gene_df[gene_df.columns].applymap(
            lambda cell: [str(item) for item in cell]
        )
        gene_df[gene_df.columns] = gene_df[gene_df.columns].applymap(
            lambda cell: ["na"] if all("nan" == item.lower() for item in cell) else cell
        )
        gene_df[gene_df.columns] = gene_df[gene_df.columns].applymap(
            lambda cell: ["na" if "nan" == item.lower() else item for item in cell]
        )

        gene_df[gene_df.columns] = gene_df[gene_df.columns].applymap(
            lambda cell: " | ".join(cell)
        )
        gene_df = gene_df[gene_df.columns]

        # annotate the passed in dataframe
        sample_df = sample_df.drop(gene_col, axis=1).join(gene_df)
        return sample_df

    def add_decipher_link(self, df):
        df["DECIPHER_LINK"] = [
            """=HYPERLINK("https://decipher.sanger.ac.uk/browser#q/%s:%s-%s")"""
            % index[0:3]
            for index, fields in df.iterrows()
        ]
