import allel
import numpy as np
import pandas as pd
import os
from collections import defaultdict
from pybedtools import BedTool
from pysam import VariantFile

# original code written by Dennis Kao: https://github.com/ccmbioinfo/crg/blob/master/SVRecords/SVGrouper.py

class SVGrouper:
    def __init__(self, vcfs, report_type, ann_fields=[]):
        def list2string(col):
            return col.apply(lambda x: ", ".join(x) if isinstance(x, list) else x)

        # key - dataframe col name : value - vcf field name
        self.index_cols = {
            "variants/CHROM": "CHROM",
            "variants/POS": "POS",
            "variants/END": "END",
            "variants/SVTYPE": "SVTYPE",
        }

        # parse vcf to return bedtools intervals, annotation fields, and samples
        all_sv, ann_df, sample_list = self._parse_sv_vcfs(
            vcfs, report_type, ann_fields=ann_fields
        )
        assert len(sample_list) == len(set(sample_list)), (
            "Duplicate sample names among input vcf's detected: %s" % sample_list
        )

        # define columns for grouped svs
        columns = ["CHROM", "POS", "END", "SVTYPE", "Ensembl Gene ID", "N_SAMPLES"]
        columns.extend(sample_list)
        columns.extend(["%s_SV_DETAILS" % s for s in sample_list])
        columns.extend(["%s_GENOTYPE" % s for s in sample_list])
        if report_type == "BND":
            columns.extend(["%s_DEPTH" % s for s in sample_list])

        self.sample_list = sample_list
        # create empty dataframe to add grouped svs
        self.df = pd.DataFrame(columns=columns).set_index(
            keys=list(self.index_cols.values())
        )
        # group (merge) svs
        self._group_sv(all_sv, report_type)
        # list set typecast is to ensure that we get a list of unique ensemble identifiers
        self.df["Ensembl Gene ID"] = self.df["Ensembl Gene ID"].apply(
            lambda gene_list: list(set(gene_list.split(",")))
            if "," in gene_list
            else [gene_list]
        )
        self.bedtool = self.make_ref_bedtool()

        # if vcfs are manta for BND report, need to make a dict of BNDs with split reads (SR) and discordant read pairs (PR) supporting the alt allele
        # this is because the allel.read_vcf function does not properly parse these fields: it only extracts the SR and PR supporting the reference allele
        if report_type == "BND":
            print("Extracting split reads and discordant reads supporting breakends...")
            bnd_dict = self._make_bnd_dict(vcfs)
            # using dict, extract SR and PR for variants in self.df
            for sample in self.sample_list:
                self.df = self._add_bnd_supp_reads(bnd_dict, sample)

        # convert listed annotations into comma separated strings
        for name in self.sample_list:
            self.df["%s_SV_DETAILS" % name] = list2string(
                self.df["%s_SV_DETAILS" % name]
            )
            self.df["%s_GENOTYPE" % name] = list2string(self.df["%s_GENOTYPE" % name])
            if report_type == "BND":
                # if there are multiple BND at one site, just take the first depth (should be same for both)
                self.df["%s_DEPTH" % name] = [
                    depth[0] if isinstance(depth, list) and len(depth) >= 1 else depth
                    for depth in (self.df["%s_DEPTH" % name].values.tolist())
                ]
                self.df["%s_DEPTH" % name] = list2string(self.df["%s_DEPTH" % name])

        # append annotation fields to final df
        self.df = self.df.join(ann_df, how="left")

    def _parse_sv_vcfs(self, vcf_paths, report_type, ann_fields=[]):
        """
        Merge all SV interval data from multiple vcf's in to a single BedTool instance

        Implementation:
            Use Panda's dataframe for some easy preprocessing, then create a BedTool from a tuple containing each row
        """

        def split_Ensembl_ids(id_list):
            new_list = []
            for id in id_list:
                if "-" in id:
                    new_list.extend(id.split("-"))
                elif "&" in id:
                    new_list.extend(id.split("&"))
                else:
                    new_list.append(id)
            return new_list

        intervals = []
        sample_names = []
        ann_dfs = []

        # CHR POS STOP needs to be first 3 columns for creation of BedTool instance
        index_fields = list(self.index_cols.keys())

        # get sample-specific ALT, GT fields, as well as breakend depth (BND_DEPTH) if this a Manta BND vcf
        # these columns will be used in the bedtools intersection
        sample_sv_fields = (
            index_fields
            + ["variants/ALT"]
            + ["calldata/GT", "variants/ANN_Gene_ID", "samples"]
        )

        if report_type == "BND":
            sample_sv_fields = sample_sv_fields + [
                "variants/BND_DEPTH",
            ]
        parse_fields = list(set(sample_sv_fields + ann_fields))

        # use read_vcf because genotype field is not picked up with vcf_to_dataframe
        for vcf_path in vcf_paths:
            vcf_dict = allel.read_vcf(
                vcf_path,
                ["*"],
                numbers={"ANN": 1000},
                transformers=allel.ANNTransformer(),
            )

            assert (
                len(vcf_dict["samples"]) == 1
            ), "%s contains 0 or more than 1 sample: %s" % (
                vcf_path,
                str(vcf_dict["samples"]),
            )
            name = vcf_dict.pop("samples")[0]
            sample_names.append(name)

            # if 'chr' in CHROM field, remove
            vcf_dict["variants/CHROM"] = [
                chrom.strip("chr") for chrom in vcf_dict["variants/CHROM"]
            ]

            # grab alt allele
            vcf_dict["variants/ALT"] = [alt[0] for alt in vcf_dict["variants/ALT"]]

            # drop un-needed fields from vcf, cannot pass in parse_fields to read_vcf() because ANN_gene_id is unknown until ANNTransformer runs
            for key in list(vcf_dict.keys()):
                if key not in parse_fields:
                    vcf_dict.pop(key)

            # remove empty strings, split on delimited characters, then join using comma
            # by default, specifying numbers=1000 creates 1000 elements, with most being empty
            vcf_dict["variants/ANN_Gene_ID"] = [
                list(filter(None, ann)) for ann in vcf_dict["variants/ANN_Gene_ID"]
            ]
            vcf_dict["variants/ANN_Gene_ID"] = [
                split_Ensembl_ids(id_list)
                if any("&" in id for id in id_list) or any("-" in id for id in id_list)
                else id_list
                for id_list in vcf_dict["variants/ANN_Gene_ID"]
            ]
            vcf_dict["variants/ANN_Gene_ID"] = [
                ",".join(list(set(id_list))) if isinstance(id_list, list) else id_list
                for id_list in vcf_dict["variants/ANN_Gene_ID"]
            ]

            # unlist Manta BND_DEPTH
            if report_type == "BND":
                vcf_dict["variants/BND_DEPTH"] = [
                    depth for depth in vcf_dict["variants/BND_DEPTH"]
                ]

            # get genotypes
            vcf_dict["calldata/GT"] = np.array(
                [
                    "HET" if 0 in gt and 1 in gt else "HOM"
                    for gt in vcf_dict.pop("calldata/GT")
                ]
            )

            df = pd.DataFrame(vcf_dict)
            df["samples"] = name

            # workaround for START > END so BedTool doesn't freak out using MalformedBedError
            # START > END is the case for TRV, INV
            s = df["variants/END"] < df["variants/POS"]
            df.loc[s, ["variants/END", "variants/POS"]] = df.loc[
                s, ["variants/POS", "variants/END"]
            ].values
            df["variants/POS"] = df["variants/POS"].astype(int)
            df["variants/END"] = df["variants/END"].astype(int)
            df = df.drop_duplicates()

            # for BND, make POS=END
            bnd = df["variants/SVTYPE"] == "BND"
            df.loc[bnd, ["variants/POS"]] = df.loc[bnd, ["variants/END"]].values
            df["variants/END"] = df["variants/END"].astype(int)
            df["variants/POS"] = df["variants/POS"].astype(int)
            df = df.drop_duplicates()

            intervals.extend(df[sample_sv_fields].itertuples(index=False))
            print(df[sample_sv_fields].itertuples(index=False))

            if ann_fields:
                ann_dfs.append(df[index_fields + ann_fields])

        ann_df = (
            pd.concat(ann_dfs)
            .astype(str)
            .rename(columns=self.index_cols)
            .set_index(list(self.index_cols.values()))
            if ann_fields
            else pd.DataFrame()
        )
        # annotations for the same SV in a vcf can have slighly differing fields (ex. SVSCORE_MEAN)
        ann_df = ann_df[~ann_df.index.duplicated(keep="first")]

        return BedTool(intervals), ann_df, sample_names

    def _group_sv(self, bedtool, report_type, reciprocal_overlap=0.5):
        already_grouped_intervals = set()

        for l in bedtool.intersect(
            bedtool, wa=True, wb=True, F=reciprocal_overlap, f=reciprocal_overlap
        ):
            if report_type == "BND":
                (
                    ref_chr,
                    ref_start,
                    ref_end,
                    ref_svtype,
                    ref_alt,
                    ref_gt,
                    ref_genes,
                    ref_name,
                    ref_bnd_depth,
                    samp_chr,
                    samp_start,
                    samp_end,
                    samp_svtype,
                    samp_alt,
                    samp_gt,
                    samp_genes,
                    samp_name,
                    samp_bnd_depth,
                ) = l
            else:
                (
                    ref_chr,
                    ref_start,
                    ref_end,
                    ref_svtype,
                    ref_alt,
                    ref_gt,
                    ref_genes,
                    ref_name,
                    samp_chr,
                    samp_start,
                    samp_end,
                    samp_svtype,
                    samp_alt,
                    samp_gt,
                    samp_genes,
                    samp_name,
                ) = l

            ref_interval = (
                ref_chr,
                ref_start,
                ref_end,
                ref_svtype,
                ref_alt,
            )
            samp_interval = (
                samp_chr,
                samp_start,
                samp_end,
                samp_svtype,
                samp_alt,
                samp_gt,
                samp_name,
            )
            if report_type == "BND":
                samp_interval = (
                    samp_chr,
                    samp_start,
                    samp_end,
                    samp_svtype,
                    samp_alt,
                    samp_gt,
                    samp_name,
                    samp_bnd_depth,
                )

            if (samp_interval not in already_grouped_intervals) and (
                ref_svtype == samp_svtype
            ):
                self._add_interval(ref_interval, ref_genes, samp_interval, report_type)
                already_grouped_intervals.add(samp_interval)

        self.df.sort_index(inplace=True)

    def _add_interval(self, ref_interval, ref_genes, samp_interval, report_type):
        if report_type == "BND":
            (
                samp_chr,
                samp_start,
                samp_end,
                samp_svtype,
                samp_alt,
                samp_gt,
                samp_name,
                samp_bnd_depth,
            ) = samp_interval
        else:
            (
                samp_chr,
                samp_start,
                samp_end,
                samp_svtype,
                samp_alt,
                samp_gt,
                samp_name,
            ) = samp_interval

        # Get reference to row
        if ref_interval[:-1] not in self.df.index:
            # make new row
            self.df.loc[ref_interval[:-1], :] = np.nan
            row = self.df.loc[ref_interval[:-1], :]
            row["N_SAMPLES"] = 0
            row["Ensembl Gene ID"] = ref_genes

            for name in self.sample_list:
                row[name] = 0
                row["%s_SV_DETAILS" % name] = []
                row["%s_GENOTYPE" % name] = []
                if report_type == "BND":
                    row["%s_DEPTH" % name] = []
        else:
            row = self.df.loc[ref_interval[:-1], :]

        # Set values for row
        if row[samp_name] == 0:
            row["N_SAMPLES"] += 1
            row[samp_name] = 1

        if samp_svtype == "BND":
            row["%s_SV_DETAILS" % samp_name].append(
                "{}:{}-{}:{}:{}".format(
                    samp_chr, samp_start, samp_end, samp_svtype, samp_alt
                )
            )
        else:
            row["%s_SV_DETAILS" % samp_name].append(
                "{}:{}-{}:{}".format(samp_chr, samp_start, samp_end, samp_svtype)
            )
        row["%s_GENOTYPE" % samp_name].append(samp_gt)
        if report_type == "BND":
            row["%s_DEPTH" % samp_name].append(samp_bnd_depth)

    def write(self, outfile_name):
        self.df.to_csv(outfile_name, sep="\t", encoding="utf-8", na_rep=".")

    def make_ref_bedtool(self):
        return BedTool(list(self.df.index.values))

    def _make_bnd_dict(self, vcfs):
        # makes a dictionary of dictionaries, e.g. {'control_DMD': {'X:40352732-40352732:BND:CCTGC[15:74508951[': {'SR': 30, 'PR': 4}}
        bnd_dict = defaultdict(dict)
        for vcf in vcfs:
            print(vcf)
            bnd_vcf = VariantFile(vcf)
            for record in bnd_vcf.fetch():
                for sample in record.samples:
                    sr = record.samples[sample]["SR"][1]
                    pr = record.samples[sample]["PR"][1]
                    name = record.samples[sample].name
                bnd = (
                    record.chrom
                    + ":"
                    + str(record.pos)
                    + "-"
                    + str(record.pos)
                    + ":"
                    + record.info["SVTYPE"]
                    + ":"
                    + "".join(
                        record.alts
                    )  # alts is a tuple because the program accounts for possibility of multiple vcfs, but we always have one-sample vcfs
                )
                bnd_dict[name][bnd] = {"SR": sr, "PR": pr}
        return bnd_dict

    def _add_bnd_supp_reads(self, bnd_dict, sample):
        df = self.df
        sr_all = []
        pr_all = []
        for index, row in df.iterrows():
            sv_details = row[f"{sample}_SV_DETAILS"]
            if sv_details == "":
                # BND was not called in sample
                sr_all.append(".")
                pr_all.append(".")
            else:
                # sv_details = sv_details.split(",")
                sr_max = 0
                pr_max = 0
                for sv in sv_details:
                    sv = sv.strip(" ")
                    sr = bnd_dict[sample][sv]["SR"]
                    pr = bnd_dict[sample][sv]["PR"]
                    if sr > sr_max:
                        sr_max = sr
                    if pr > pr_max:
                        pr_max = pr
                sr_all.append(sr_max)
                pr_all.append(pr_max)

        df[f"{sample}_SR"] = sr_all
        df[f"{sample}_PR"] = pr_all
        return df
