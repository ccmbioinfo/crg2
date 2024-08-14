import sys, os, re
from collections import namedtuple

'''
merge output files of EHDN after running through DBSCAN script into one file
get values for outlier, repeat size, a1000g_freq from EHdn.expansion.<date>.tsv 
and merge with merged.rare.expansion.<date>.tsv for ANNOVAR annotation
manifest_1000G.txt: 
    EHDN manifest files (three-column:sampleid\tcase/control\tpath to EHDN JSON)
    you can also provide a single column file with 1000G sample ids; this is uesd
    to remove rows with only 1000G samples as outlier
python format_for_annovar.py EHdn.expansions.2021-02-18.tsv merged.rare.expansions.2021-02-18.tsv manifest_1000G.txt
'''

#not using pandas.merge because the two files have some duplicate lines (chr#start#end)
#merged.rare.expansion has altered intervals (merging nearby motifs) so won't exactly match those
#in EHdn.expansion.<dat>.tsv file entries

ehdn_out = sys.argv[1]
dbscan_out = sys.argv[2]
g1k_manifest = sys.argv[3]
date = ehdn_out.split(".")[-2]
d = os.path.dirname(dbscan_out)
outfile = os.path.join(d,"merged.rare.EHdn.expansion." + date + ".tsv")
samples = []
print("input files: ", ehdn_out, dbscan_out, g1k_manifest)
print("out:", outfile)

#read 1000G sample names from manifest file
g1k = []
with open(g1k_manifest) as f:
    for i in f:
        i = i.strip("\n")
        if re.split("\t ' '",i)[0]:
            g1k.append(re.split("\t|' '",i)[0])
        else:
            g1k.append(i)
all_g1k_sample = lambda x: set(x).difference(g1k) 

#read DBSCAN output
dbscan = namedtuple("dbscan", "chr start end motif outliers key")
merged_exp = []
with open(dbscan_out) as f:
    f.readline()
    for i in f:
        i = i.strip().split("\t")
        key = "#".join(i[0:3])
        outlier = re.split(";",i[4])
        if all_g1k_sample(outlier):
            merged_exp.append(dbscan(chr=i[0],start=i[1],end=i[2],motif=i[3],outliers=i[4],key=key))
            for o in outlier:
                if not o in g1k:
                    samples.append(o)

samples = list(set(samples))
print("samples = ", samples)
#read EHDN output
ehdn = namedtuple("ehdn", "motif outliers size ref chr start end a1000g key")
ehdn_exp = []
with open(ehdn_out) as f:
    f.readline()
    for i in f:
        i = i.strip().split("\t")
        key = "#".join(i[5:8])
        ehdn_exp.append(ehdn(motif=i[1],outliers=i[2],size=i[3],ref=i[4],chr=i[5],start=i[6],end=i[7],a1000g=i[8],key=key))
ehdn_exp        

#match reacords from DBSCAN and EHDN to retrieve all info
#DBSCAN repeats are merged based on locations, so may not match original EHDN output

def split_outlier_size(outlier,size,samples):
    outlier_l = re.split(";",outlier) 
    size_l = re.split(";",size)
    if len(size_l) == len(outlier_l):
        len_samples = len(samples)
        outlier_samples = ["0"]*(len_samples+1)
        outlier_sizes = ["0"]*len_samples
        for i,s in enumerate(samples):
            if s in outlier_l:
                ind = outlier_l.index(s)
                outlier_samples[i] = "1"
                outlier_sizes[i] = size_l[ind]
            else:
                if s in g1k:
                    outlier_samples[len_samples] += 1
    return outlier_samples, outlier_sizes

found, allkey, out, seen = [], [], [], []
ref, alt = "0", "-"
for i in merged_exp:
    allkey.append(i)
    for j in ehdn_exp:
        if i.key == j.key:
            if i.outliers == j.outliers:
                found.append(i)
                outliers, sizes = split_outlier_size(i.outliers, j.size, samples)
                print(outliers, sizes, i.key)
                # out.append([i.chr, i.start, i.end, ref, alt, i.motif, i.outliers, i.key, j.size, j.a1000g])
                # out.append([i.chr, i.start, i.end, ref, alt, i.motif, outliers, sizes, i.key, j.a1000g])
                o = [i.chr, i.start, i.end, ref, alt, i.motif ] +  outliers + sizes + [i.key, j.a1000g]
                out.append(o)
                print(o)
            # seen.append(i)
                continue


notfound = list(set(allkey).difference(found))
within = lambda a,b,c,d: ( c <= a and a <= d ) or ( c <= b and b <= d)
for i in notfound:
    for j in ehdn_exp:
        if i.chr == j.chr and i.motif == j.motif and i.outliers == j.outliers:
            r1 = range(int(i.start),int(i.end))
            r2 = range(int(j.start),int(j.end))
#            if i.start == j.start or i.end == j.end or len(set(r1).intersection(r2)) > 5:
            if within(i.start,i.end,j.start,j.end) or within(j.start,j.end,i.start,i.end):
                # out.append([i.chr, i.start, i.end, ref, alt, i.motif, i.outliers, i.key, j.size, j.a1000g])
                outliers, sizes = split_outlier_size(i.outliers, j.size, samples)
                o = [i.chr, i.start, i.end, ref, alt, i.motif ] +  outliers + sizes + [i.key, j.a1000g]
                # out.append([i.chr, i.start, i.end, ref, alt, i.motif, outliers, sizes, i.key, j.a1000g])
                out.append(o)
                print(o)
#               seen.append(i)
                continue

columns = ["chr", "start", "end", "ref", "alt", "motif", ]
columns += ["outlier." + s for s in samples]
columns += ["1000G_outlier_count"]
columns += ["size." + s for s in samples]
columns += ["key", "1000G_freq_perc"]
with open(outfile, "w")  as f:
    for i in columns:
        f.writelines("%s\t"%i)
    f.writelines("\n")
    for i in out:
        for j in i:
            f.writelines("%s\t"%j)
        f.writelines("\n")
