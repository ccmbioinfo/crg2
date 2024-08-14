'''
reads in disease threshold annotated EH results
splits annotation into different columns
fills outlier column with sample name if GT is in disease threshold
adds mean, std, median GT sizes from 1000Genomes EH runs (static file)
writes o/p as excel file
'''

import sys, os, re
from collections import namedtuple
from xlsxwriter.workbook import Workbook

annot_tsv = sys.argv[1] #output from add_gene+threshold_to_EH_column_headings2.py
xlsx = sys.argv[2] #output xlsx filename

def outlier_gt(threshold, gt_dict):
    outlier = []
    for sample in gt_dict:
        threshold = re.sub(u'\u2014','-',str(threshold))
        #threshold = str(threshold.strip())
        x = gt_dict[sample]
        y = x.split("/") if "/" in x else x
        if not "." in y:
            y = max(list(map(int,y))) if isinstance(y,list) else int(y)
        else:
            return " "
        if y:
            if "(" in threshold or " " in threshold:
                continue
            elif ">" in threshold:
                if "-" in threshold:
                    condition = int(threshold.split(">")[1].split("-")[0])
                else:
                    condition = int(threshold.split(">")[1])
                    if y >= condition:
                        outlier.append(sample)
            elif "-" in threshold:
                condition = list(map(int, threshold.split("-")))
                if y in range(condition[0], condition[1]):
                    outlier.append(sample)
            elif "," in threshold:
                condition = list(map(int, threshold.split(",")))
                if y in condition:
                    outlier.append(sample)
            else:
                if threshold.isdigit() and y >= int(threshold):
                    outlier.append(sample)    
    return ",".join(outlier) if outlier else ' '
    


EH = namedtuple('EH', 'pos, motif, gene, size, gt, mean, std, median')

G1K = {}
with open(os.path.expanduser("~/crg/1000G_EH_v1.0.tsv")) as f:
    for i in f:
        if not i.startswith("#"):
            annot, mean, std, median = i.strip().split("\t")
            if not annot in G1K:
                G1K[annot] = [mean, std, median]

eh_gt = {}
with open(annot_tsv) as f:
    annot = {e:l for e,l in enumerate(f.readline().strip().split("\t")[1:]) }
    for i in f:
        i = i.strip().split("\t")
        sample = i[0]
        if not sample in eh_gt:
            eh_gt[sample] = []

        for  e,l in enumerate(i[1:]):
            pos, motif, gene, size = annot[e].rsplit(":",3)
            gt = l
            if annot[e] in G1K:
                mean, std, median = G1K[annot[e]]
            else:
                mean, std, median = "NA", "NA", "NA"
            eh_gt[sample].append(EH(pos,motif,gene,size,gt,mean,std,median))

trf = {}
for i in eh_gt:
    for items in eh_gt[i]:
        if not items.pos in trf:
            trf[items.pos] = {}
        trf[items.pos].update({i:items})

samples = list(eh_gt.keys())

xlsx = xlsx if ".xlsx" in xlsx else xlsx + ".xlsx"
workbook = Workbook(xlsx)

worksheet = workbook.add_worksheet("Allsamples")
header = ["#location", "repeat motif", "gene", "disease threshold"] + ["GT."+i for i in samples] + ["1000G_mean", "1000G_std", "1000G_median", "outlier"]
worksheet.write_row(0, 0, header)
row = 1
for i in trf:
        info = trf[i][samples[0]]
        content = [info.pos, info.motif, info.gene, info.size]
        content += [ trf[i][s].gt for s in samples ]          
        content += [info.mean, info.std, info.median] 
        gt = { s:trf[i][s].gt for s in samples }
        outlier = outlier_gt(info.size, gt)
        content += [ outlier ]

        worksheet.write_row(row, 0, content)
        row += 1
workbook.close()

    


        
