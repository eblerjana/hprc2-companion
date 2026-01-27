import sys, argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import upsetplot
import pandas as pd
import gzip

outname = sys.argv[1]

all_rare_svs_covered = []
all_svs_covered = []
all_sv_errors = []

all_small_covered = []
all_small_errors = []


for filename in sys.stdin:
	df = pd.read_csv(filename.strip(), sep='\t')
	columnnames = list(df)
	print(columnnames)

	polished = columnnames[9]
	unpolished = columnnames[10]

	assert '-polished' in polished
	assert 'unpolished' in unpolished

	all_rare_svs_covered.append( len(df[ ~(df["overlaps_SV"]) & (df["in_assembly"]) & (df[polished]) & (df["variant_length"] >= 50) ]) - len(df[ ~(df["overlaps_SV"]) & (df["in_assembly"]) & (df[unpolished]) & (df["variant_length"] >= 50) ]) )
	all_svs_covered.append( len(df[ (df["in_assembly"]) & (df[polished]) & (df["variant_length"] >= 50) ]) - len(df[ (df["in_assembly"]) & (df[unpolished]) & (df["variant_length"] >= 50) ]) )
	all_sv_errors.append( len(df[ ~(df["in_assembly"]) & (df[polished]) & (df["variant_length"] >= 50) ]) - len(df[ ~(df["in_assembly"]) & (df[unpolished]) & (df["variant_length"] >= 50) ]) )
	all_small_covered.append( len(df[ (df["in_assembly"]) & (df[polished]) & (df["variant_length"] < 50) ]) - len(df[ (df["in_assembly"]) & (df[unpolished]) & (df["variant_length"] < 50) ]) )
	all_small_errors.append( len(df[ ~(df["in_assembly"]) & (df[polished]) & (df["variant_length"] < 50) ]) - len(df[ ~(df["in_assembly"]) & (df[unpolished]) & (df["variant_length"] < 50) ]) )


print(all_rare_svs_covered)
print(all_svs_covered)
print(all_sv_errors)

print(all_small_covered)
print(all_small_errors)


plt.style.use('tableau-colorblind10')
plt.rcParams["font.family"] = "Nimbus Sans"

# create upset plot
with PdfPages(outname) as outpdf:

	values = [all_rare_svs_covered, all_svs_covered, all_sv_errors, all_small_covered, all_small_errors]
	labels = ['rare SVs\ncovered', 'total\nSVs covered', 'SV errors', 'total\nSNPs+indels covered', 'SNP+indel\nerrors']

	fig, ax = plt.subplots()
	plt.ylabel('Variant Count Change After Polishing')
	plt.boxplot(values, tick_labels = labels)
#	ax.yaxis.grid(True)
#	ax.xaxis.grid(True)
	plt.axhline(y=0, color='black', linewidth=0.2)
	ax.set_yscale('symlog')
	ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

	plt.tight_layout()
	outpdf.savefig()
	plt.close()


