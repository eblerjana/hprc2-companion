import sys
import argparse
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

def parse_afs(filename):
	afs = []
	for line in open(filename, 'r'):
		if line.startswith('#chromosome'):
			continue
		fields = line.strip().split()
		if fields[8] == "None":
			continue
		for af in fields[8].split(','):
			afs.append(float(af))
	return afs

def plot_missed_afs(files, outfile, samples, sampling_sizes):
	fig, axes = plt.subplots(figsize=(35,5), nrows=1, ncols=len(sampling_sizes))
	for i,s in enumerate(sampling_sizes):
		# create data frame containing data for all samples
		afs = []
		for sample in samples:
			afs.append(parse_afs(files[(sample, s)]))

		if len(sampling_sizes) > 1:
			axes[i].violinplot(afs)
			axes[i].set_xticks([a+1 for a in range(len(samples))])
			axes[i].set_xticklabels(samples)
			axes[i].set_title('sampling parameters: ' + s)
			axes[i].set_ylabel('allele frequencies of missed alleles')
		else:
			axes.violinplot(afs)
			axes.set_xticks([a+1 for a in range(len(samples))])
			axes.set_xticklabels(samples)
			axes.set_title('sampling parameters: ' + s)
			axes.set_ylabel('allele frequencies of missed alleles')


	fig.tight_layout()
	plt.savefig(outfile)



if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='plot_missed_afs.py', description = __doc__)
	parser.add_argument('-tsv', metavar='TSVS', nargs='+', required=True, help='Panel stats computed by evaluate_panel.py script.')
	parser.add_argument('-outfile', metavar='OUTFILE', required=True, help='Name of the output file (PDF format).')
	args = parser.parse_args()

	sampling_sizes = []
	files = {}
	samples = set([])
	sampling_sizes = set([])
	for t in args.tsv:
		sampling_size = t.split('/')[-1].split('.tsv')[0].split('_')[-1]
		sample = t.split('/')[-1].split('.tsv')[0].split('_')[-2]
		files[(sample, sampling_size)] = t
		samples.add(sample)
		sampling_sizes.add(sampling_size)
	samples = sorted(list(samples))
	sampling_sizes = sorted(list(sampling_sizes), key=lambda s: int(s.split('-')[0]))


	plot_missed_afs(files, args.outfile, samples, sampling_sizes)
