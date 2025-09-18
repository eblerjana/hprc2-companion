import sys, argparse
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import upsetplot
import pandas as pd
import gzip

parser = argparse.ArgumentParser(prog='plot_upset.py', description=__doc__)
parser.add_argument('-c', '--config', required=True, help='Config file: assemblies.tsv used for PAV.')
parser.add_argument('-o', '--outname', required=True, help='output name prefix.')
parser.add_argument('-t', '--threshold', required=True, type=float, help='Threshold to consider for plotting SVs.')
parser.add_argument('-v', '--variants', required=True, type=string, help='File with TP variant ids.')
args = parser.parse_args()


callsetname_to_index = defaultdict(lambda: [])
match_ids = set([])

# parse IDs of TP variant matches
for line in gzip.open(args.v, 'rt'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	match_ids.add(fields[2])

# parse PAV assembly.tsv file to determine order of callsets
for line in open(args.config, 'r'):
	if line.startswith('NAME'):
		continue
	fields = line.strip().split()
	for i,field in enumerate(fields[1:]):
		name = field.split('_')[0]
		callsetname_to_index[name].append(i)

callsetnames = sorted(callsetname_to_index.keys())

with open(args.outname + '.tsv', 'w') as outtsv:

	# print header
	header_line = ["ID", "variant_type", "variant_length", "overlaps_bubble", "overlaps_SV"] + ['in_' + c for c in callsetnames] + ['GT_' + c for c in callsetnames]
	outtsv.write("\t".join(header_line) + "\n")

	# print variant lines
	for line in sys.stdin:
		if line.startswith('#'):
			continue
		fields = line.strip().split()
		if fields[6] != "PASS":
			continue
		if "X" in fields[0] or "Y" in fields[0]:
			continue
		info_field = {f.split('=')[0] : f.split('=')[1] for f in fields[7].split(';') if '=' in f}
		assert "SVTYPE" in info_field
		assert "ID" in info_field
		vartype = info_field["SVTYPE"]
		varid = fields[2]
		variant_length = 1 if "SNV" in varid else int(varid.split('-')[-1].split('.')[0])
		gt_index = fields[8].split(':').index("GT")
		alleles = fields[9].split(':')[gt_index].split('|')
		overlaps_regions = float(fields[10]) > 0.0

		# print line
		overlaps_sv = varid in match_ids 
		line_to_print = [varid, vartype, str(variant_length), str(overlaps_regions), str(overlaps_sv)]

		for callset in callsetnames:
			assert len(callsetname_to_index[callset]) == 2
			allele1 = alleles[callsetname_to_index[callset][0]]
			allele2 = alleles[callsetname_to_index[callset][1]]
			if (allele1 == "1") or (allele2 == "1"):
				line_to_print.append("True")
			else:
				line_to_print.append("False")

		for callset in callsetnames:
			assert len(callsetname_to_index[callset]) == 2
			allele1 = alleles[callsetname_to_index[callset][0]]
			allele2 = alleles[callsetname_to_index[callset][1]]
			line_to_print.append('|'.join([allele1, allele2]))
	
		outtsv.write("\t".join(line_to_print) + "\n")


# create upset plot
with PdfPages(args.outname + ".pdf") as outpdf:
	df = pd.read_csv(args.outname + '.tsv', sep='\t')
	columnnames = ['in_' + c for c in callsetnames]

	combinations = [columnnames]
	for callset in callsetnames:
		if "assembly" in callset:
			continue
		combinations.append(['in_assembly', 'in_' + callset])

	plt.figure()

	for combination in combinations:

		# plot overlaps for all variants
		fig = plt.figure()
		variants = df.groupby(by=combination).size()
		upsetplot.plot(variants, sort_by='cardinality', show_counts='%d')
		fig.tight_layout()
		outpdf.savefig()
		plt.close()
	
		# plot overlaps for small variants
		fig = plt.figure()
		df_sub = df[df["variant_length"] < args.threshold]
		variants = df_sub.groupby(by=combination).size()
		upsetplot.plot(variants, sort_by='cardinality', show_counts='%d')
		plt.suptitle('small variants (< ' + str(args.threshold) +' bp)')
		fig.tight_layout()
		outpdf.savefig()
		plt.close()

		# plot overlaps for SVs 
		fig = plt.figure()
		df_sub = df[df["variant_length"] >= args.threshold]
		variants = df_sub.groupby(by=combination + ['overlaps_SV']).size()
		upsetplot.plot(variants, sort_by='cardinality', show_counts='%d')
		plt.suptitle('SVs (>= ' + str(args.threshold) + 'bp)')
		fig.tight_layout()
		outpdf.savefig()
		plt.close()


	# print some statistics
	
	for callset in callsetnames:
		if "assembly" in callset:
			continue
		print("-------------------------------------")
		print("SV statistics for " + callset + ":")

		rare_count = len(df[ ~(df["overlaps_bubble"]) & (df["in_assembly"]) & (df["in_" + callset]) & (df["variant_length"] >= args.threshold) ])
		total = len(df[ ~(df["overlaps_bubble"]) & (df["in_assembly"]) & (df["variant_length"] >= args.threshold) ])
		print("rare ground truth SVs covered: " + str(rare_count))
		print("rare ground truth SVs: " + str(total))
		print("")

		truth_count = len(df[ (df["in_assembly"]) & (df["in_" + callset]) & (df["variant_length"] >= args.threshold) ])
		total = len(df[ (df["in_assembly"]) & (df["variant_length"] >= args.threshold) ])
		print("ground truth SVs covered: " + str(truth_count))
		print("ground truth SVs: " + str(total))
		print("")

		error_count = len(df[ ~(df["in_assembly"]) & (df["in_" + callset]) & (df["variant_length"] >= args.threshold) ])
		total = len(df[ (df["in_" + callset]) & (df["variant_length"] >= args.threshold) ])
		print("error SVs in callset: " + str(error_count))
		print("total SVs in callset: " + str(total))
		print("")

		print("small variant statistics for " + callset + ":")

		truth_count = len(df[ (df["in_assembly"]) & (df["in_" + callset]) & (df["variant_length"] < args.threshold) ])
		total = len(df[ (df["in_assembly"]) & (df["variant_length"] < args.threshold) ])
		print("ground truth variants covered: " + str(truth_count))
		print("ground truth variants: " + str(total))
		print("")

		error_count = len(df[ ~(df["in_assembly"]) & (df["in_" + callset]) & (df["variant_length"] < args.threshold) ])
		total = len(df[ (df["in_" + callset]) & (df["variant_length"] < args.threshold) ])
		print("error variants in callset: " + str(error_count))
		print("total variants in callset: " + str(total))
		print("")
		print("-------------------------------------")
