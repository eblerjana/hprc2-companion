import sys
import argparse
import pandas as pd
from collections import defaultdict
	
def print_summary(tsv, sample, sampling_size):
	df = pd.read_csv(tsv, sep='\t')
	chromosomes = ['whole_genome'] + list(df['#chromosome'].unique())
	header = list(df.columns)
	hap_columns = [c for c in header if c.startswith('Recombination_')]
	nr_paths = len(hap_columns)
	print('\t'.join(['sample', 'sampling_size', 'chromosome'] + hap_columns))
	df = df[['#chromosome'] + hap_columns]
	
	for chrom in chromosomes:
		counts = [0] * nr_paths
		df_region = df if chrom == 'whole_genome' else df[df['#chromosome'] == chrom]
		for i,hap in enumerate(hap_columns):
			recomb_count = df_region[hap].value_counts().get(1, 0)
			counts[i] = recomb_count
		print('\t'.join([sample, sampling_size, chrom] + [str(c) for c in counts]))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='evaluate_recombination.py', description=__doc__)
	parser.add_argument('-tsv', metavar='TSV', help='TSV with recombination stats.', required=True)
	parser.add_argument('-sample', metavar='SAMPLE', help='sample name.', required=True)
	parser.add_argument('-size', metavar='SIZE', help='sampling size', required=True)
	args = parser.parse_args()

	print_summary(args.tsv, args.sample, args.size)
