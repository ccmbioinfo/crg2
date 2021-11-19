import yaml, os, argparse, re
import numpy as np
from xlsxwriter.workbook import Workbook

'''
Replace existing 3 scripts with single script, by writing new functions here
1. Parse EH JSON output
2. Extract location, repeat unit and genotype of repeat
3. Annotate with gene name, disease threshold
4. Automatically check if GT exceeds the disease threshold or 1000G mean+3sd
	4.1 if yes, create a outlier column with sample names for that location
5. Write output to xlsx format
'''

def read_eh_json(filename, locations):

	'''
	Read JSON from Expansion Hunter
	Extract genotype info
	'''

	inp = yaml.safe_load(open(filename))
	eh = {}
	#EH profile
	eh_results = inp['LocusResults']
	catalog = eh_results.keys()	      
	for i in catalog:           
		#for EH v3.0 -> 5.0 dict keys are structured this way; inspect JSON output and change if you get any key errors
		loc, repeat, gt = [ eh_results[i]['Variants'][i][k] for k in ['ReferenceRegion', 'RepeatUnit', 'Genotype'] ]
		#remove chr prefix if present
		if loc.startswith("chr"):
			loc = loc.replace("chr","")
		key = loc + ":" + repeat
		if not key in eh:
			eh[key] = gt
		if not key in locations:
			locations[key] = 1
	return eh

def parse_disease_threshold(threshold, x):

	'''
	Parse disease threshold text column to integer comparison with genotype
	return True when genotype >= disease_threshold
	'''
	condition = ""
	if "(" in threshold or " " in threshold:
		return False
	elif ">" in threshold:
		if "-" in threshold:
			# condition = int(threshold.split(">")[1].split("-")[0])
			condition = min(list(map(int, threshold.split(">")[1].split("-"))))
			if x >= condition:
				return True
		else:
			condition = int(threshold.split(">")[1])
			if x >= condition:
				return True
	elif "-" in threshold:
		condition = min(list(map(int, threshold.split("-"))))
		if x >= condition:
			return True

	elif "," in threshold:
		condition = list(map(int, threshold.split(",")))
		if x in condition:
			return True
	else:
		if threshold.isdigit() and x >= int(threshold):
			return True
	return False

def check_outlier(disease_threshold, sample_gt, g1k=None):

	'''
	return True if genotype >= disease_threshold or >= 1000G(mean+3sd)
	keeping it as "OR" condition otherwise
	since EH does not accurately genotype FMR1, gennotyped size for it 
	will never exceed the disease threshold but will be marked as
	outlier by 1000G(mean+3sd). Outlier feature is manily used for
	screening for further investigation.
	'''

	x = sample_gt.split("/") if "/" in sample_gt else sample_gt

	if not "." in x:
		x = max(list(map(int,x))) if isinstance(x,list) else int(x)
	if x:
		if parse_disease_threshold(disease_threshold, x):
			return True

		else:
			if g1k and x >= ( g1k["mean"] + 3 * g1k["std"]):
				return True
			else:
				return False
	return False


def write_xlsx(xlsx, annot_str, repeat_annot, g1k=None):

	workbook = Workbook(xlsx)
	worksheet = workbook.add_worksheet("Allsamples")
	header = ["#location", "repeat motif", "gene", "disease threshold"] + ["GT."+ i.split(".")[0] for sample in annot_str] 
	if g1k:
		header += ["1000G_mean", "1000G_std", "1000G_median", "outlier"]
	else:
		header += ["outlier"]
	worksheet.write_row(0, 0, header)
	row = 1
	for l in locations:
		chrom, start, end, motif = re.split(":|-",l)
		gene, disease_threshold, motif = repeat_annot[l]
		column1 = chrom + ":" + start + "-" + end
		content = [column1,motif,gene,disease_threshold]
		outliers = []
		for sample in annot_str:
			content += [annot_str[sample][l]["gt"]]
			if annot_str[sample][l]["outlier"]:
				outliers.append(sample.split(".")[0])
		outliers = ",".join(outliers) if outliers else ""
		if not g1k or not l in g1k:
			content += [outliers]
		else:
			content += [outliers, g1k[l]["mean"], g1k[l]["std"], g1k[l]["median"]]
		worksheet.write_row(row, 0, content)
		row += 1
	workbook.close()	


if __name__ == '__main__':

	description = "Parse JSON output(s) from Expansion Hunter, annotate, detect outlier, and output as Excel report"
	parser = argparse.ArgumentParser(description=description)
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument("-j", "--json", nargs="+", help="One or more JSON output from EH")
	group.add_argument("-f", "--file", help="Single column file listing the absolute path of EH JSON output")
	parser.add_argument("-a", "--annot", required=True, help="TSV file with pathogenic repeat annotation")
	parser.add_argument("-k", "--k1g", required=False, help="Three-column TSV with stats from 1000G on the same catalog above")
	parser.add_argument("-g", "--gen", default="hg19", help="set hg19 or hg38 coordinates to use")
	parser.add_argument("-o", "--out", required=True, help="output Excel filename")
	args = parser.parse_args() 
	
	locations = {} #catalog locations
	
	#list of json files on command-line
	if args.json:
		str = {os.path.basename(i):read_eh_json(i,locations) for i in args.json }
	
	#file with list of json file paths
	elif args.file:
		str = {}
		with open(args.file) as f:
			for i in f:
				i = i.strip("\n")
				if os.path.isfile(i):
					sample = os.path.basename(i)
					if not sample in str:
						str[sample] = read_eh_json(i, locations)
				else:
					print("File {} not found. Exiting!".format(i))
					exit
	else:
		print("Inputs must be passed with either -j <space-sep JSON files> or -f <file with JSON file paths>. Exiting!")
		exit

	#read annotation TSV (file passed to us by Brett Trost)
	with open(args.annot) as f:
		repeat_annot = {}
		for i in f:
			if not i.lower().startswith("reference"):
				i = i.strip().split("\t")
				if len(i) < 10:
					pmid, category, disorder, gene, threshold, motif, hg19, hg38, mask = i
				else:
					pmid, category, disorder, gene, threshold, motif, hg19, hg38, mask, note = i					
				gen_version = args.gen if args.gen else "hg19"
				coord = hg38.replace("chr","") if gen_version == "hg38" else hg19.replace("chr","")
				if "," in motif:
					for m in motif.split(","):
						key = coord + ":" + m
						if not key in repeat_annot and mask != "TRUE":
							repeat_annot[key] = [gene, threshold, m]
				else:
					key = coord + ":" + motif
					if not key in repeat_annot and mask != "TRUE":
						repeat_annot[key] = [gene, threshold, motif]
		print("repeat annotation unmasked read: {}".format(len(repeat_annot)))

	#annotate existing str dict
	# inter = list(set(locations.keys()).intersection(repeat_annot.keys()))
	# print("overlapping annotation and json locations = {}".format(len(inter)))


	#Find outlier based on disease threshold or 1000G stats
	g1k = {}
	if args.k1g:
		k1g = np.genfromtxt(args.k1g, delimiter="\t", dtype='unicode')
		g1k = { i[0].rsplit(":",2)[0]:{"mean": float(i[1]), "median": float(i[2]), "std": float(i[3])}  for i in k1g}
		
	annot_str = {}
	for i in locations:
		if i in repeat_annot:
			gene, threshold, motif = repeat_annot[i]
			for sample in str:
				gt = str[sample][i]
				if i in g1k:
					if gene == "FXN":
						print("FXN = ", gene, threshold, motif,gt, sample, check_outlier(threshold, gt, g1k[i] ))
					outlier = check_outlier(threshold, gt, g1k[i])
				else:
					outlier = check_outlier(threshold, gt)
				if not sample in annot_str:
					annot_str[sample] = {}
				annot_str[sample][i] = {"gt": gt, "outlier": outlier}
	
	xlsx = args.out
	if not xlsx.endswith(".xlsx"):
		xlsx += ".xlsx"

	if g1k:
		write_xlsx(xlsx, annot_str, repeat_annot, g1k)
	else:
		write_xlsx(xlsx, annot_str, repeat_annot, g1k)