import sys
import argparse
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from collections import defaultdict


def plot_bubble_stats(filename, outfile):
	df = pd.read_csv(filename, sep='\t')
	chromosomes = list(df['chromosome'].unique())
	print("Processing chromosomes: " + ','.join(chromosomes))

	categories = ['all_bubbles', 'bubbles>=50bp', 'bubbles<50bp', 'multiallelic_bubbles', 'biallelic_bubbles', 'UK>0', 'UK=0']
	sampling_sizes = df['sampling_size'].unique()
	samples = df['sample'].unique()
	sample_to_index = {s : i for i,s in enumerate(samples)}

	with PdfPages(outfile) as pdf:
		for chrom in chromosomes:
			fig, axes = plt.subplots(figsize=(35,18), nrows=4, ncols=len(categories))
			for i,genotype in enumerate(['all', 'absent', 'heterozygous', 'homozygous']):
				for j,category in enumerate(categories):
					to_plot = {}
					for size in sampling_sizes:
						to_plot[size] = [0.0] * len(samples)
						df_subset = df[ (df['chromosome'] == chrom) & (df['region'] == category) & (df['true_genotype'] == genotype) & (df['sampling_size'] == size) ]
						for sample, covered in zip(df_subset['sample'], df_subset['covered_alleles[%]']):
							to_plot[size][sample_to_index[sample]] = covered
					df_figure = pd.DataFrame(to_plot, index=list(samples))
					a = df_figure.plot.line(ax=axes[i,j], ylabel='ground truth alleles covered [%]',  title = genotype + ', ' + category, rot=90, style='.-', ylim=[46, 104])
			plt.setp(axes, xticks=[i for i in range(len(samples))], xticklabels=list(samples))
			fig.suptitle(chrom, size=16)
			fig.tight_layout(rect=[0, 0.03, 1, 0.95])
			pdf.savefig()
			plt.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='plot_panel_stats.py', description = __doc__)
	parser.add_argument('-tsv', metavar='TSV', required=True, help='Bubble statistics computed by evaluate_panel.py script (TSV format).')
	parser.add_argument('-outfile', metavar='OUTFILE', required=True, help='Name of the output file (PDF format).')
	args = parser.parse_args()

	plot_bubble_stats(args.tsv, args.outfile)
