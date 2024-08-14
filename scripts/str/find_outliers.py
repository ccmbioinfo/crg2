import sys, numpy as np
from matplotlib import pyplot as plt



infile = sys.argv[1]
outfile = sys.argv[2]
tr = {}
rec = []
with open(infile) as f:
    for  i in f:
        i = i.strip().split("\t")
        if any(s in i for s in ['motif', 'zscore', 'end', 'start']):
            header = i
        else:
            depth = i[-1]
            if "," in depth:
                j = depth.split(",")
                for item in j:
                    k = item.split(":")[0]
                    if not k in tr:
                        tr[k] = 0
                    tr[k] += 1
            else:
                k = depth.split(":")[0]
                if not k in tr:
                        tr[k] = 0
                tr[k] += 1

allcounts = [ counts for sample, counts in tr.items() if not sample.startswith(("HG","NA"))]
stats = lambda l: (np.mean(l), np.std(l), np.median(l))
allcounts_stats = stats(allcounts)
not_outlier = lambda l, c: (l[0] - 3 * l[1]) < c and (l[0] + 3 * l[1] ) > c
with open(outfile, "w") as f:
    for sample, counts in tr.items():
        if not sample.startswith(("HG","NA")):
            if not not_outlier(allcounts_stats, counts):
                f.writelines("%s"%sample)