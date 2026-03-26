import sys
import argparse

def parse_concordance(filename):
	concordance = None
	typed = None
	concordance_absent = None
	concordance_het = None
	concordance_hom = None
	for line in open(filename, 'r'):
		if line.startswith('weighted_concordance'):
			concordance = float(line.split()[-1]) * 100.0
		if line.startswith('typed'):
			typed = float(line.split()[-1]) * 100.0
		if line.startswith('genotype_concordance_absent'):
			concordance_absent = float(line.split()[-1]) * 100.0
		if line.startswith('genotype_concordance_het'):
			concordance_het = float(line.split()[-1]) * 100.0
		if line.startswith('genotype_concordance_hom'):
			concordance_hom = float(line.split()[-1]) * 100.0
	assert concordance is not None
	assert typed is not None
	assert concordance_absent is not None
	assert concordance_het is not None
	assert concordance_hom is not None
	return concordance, concordance_absent, concordance_het, concordance_hom, typed


def parse_precision_recall_vcfeval(filename):
	precision = None
	recall = None
	fscore = None
	for line in open(filename, 'r'):
		if 'None' in line:
			fields = line.strip().split()
			assert len(fields) == 8
			assert fields[0] == 'None'
			precision = float(fields[5]) if fields[5] != "NaN" else 0.0
			recall = float(fields[6]) if fields[6] != "NaN" else 0.0
			fscore = float(fields[7]) if fields[7] != "NaN" else 0.0
		if line.startswith('0 total baseline variants'):
			precision = 0.0
			recall = 0.0
			fscore = 0.0
	assert precision is not None
	assert recall is not None
	assert fscore is not None
	return precision, recall, fscore


def parse_precision_recall_truvari(filename):
	precision = None
	recall = None
	fscore = None
	concordance = None
	for line in open(filename, 'r'):
		if "{" in line:
			continue
		if "}" in line:
			continue
		if line == "":
			continue
		fields = line.strip().split()
		if "precision" in fields[0]:
			assert precision is None
			precision = float(fields[-1][:-1]) if fields[-1][:-1] != "null" else 0.0
		if "recall" in fields[0]:
			assert recall is None
			recall = float(fields[-1][:-1]) if fields[-1][:-1] != "null" else 0.0
		if "f1" in fields[0]:
			assert fscore is None
			fscore = float(fields[-1][:-1]) if fields[-1][:-1] != "null" else 0.0
		if "gt_concordance" in fields[0]:
			assert concordance is None
			concordance = float(fields[-1][:-1]) * 100.0
	assert precision is not None
	assert recall is not None
	assert fscore is not None
	assert concordance is not None
	return precision, recall, fscore, concordance




def collect_concordances(files, outname):
	with open(outname, 'w') as outfile:
		outfile.write('\t'.join(['sample', 'weighted_genotype_concordance', 'typed_variants', 'concordance_absent', 'concordance_het', 'concordance_hom']) + '\n')
		for f in sorted(files):
			sample = f.split('/')[-4]
			concordance, concordance_absent, concordance_het, concordance_hom, typed = parse_concordance(f)
			outfile.write('\t'.join([sample, str(concordance), str(typed), str(concordance_absent), str(concordance_het), str(concordance_hom)]) + '\n')


def collect_precision_recall(files, outname, method):
	with open(outname, 'w') as outfile:
		header = ['sample', 'precision', 'recall', 'fscore', 'gt_concordance'] if method == "truvari" else ['sample', 'precision', 'recall', 'fscore']
		outfile.write('\t'.join(header) + '\n')
		for f in sorted(files):
			sample = f.split('/')[-4]
			if method == "vcfeval":
				precision, recall, fscore = parse_precision_recall_vcfeval(f)
				outfile.write('\t'.join([sample, str(precision), str(recall), str(fscore)]) + '\n')
			else:
				precision, recall, fscore, conc = parse_precision_recall_truvari(f)
				outfile.write('\t'.join([sample, str(precision), str(recall), str(fscore), str(conc)]) + '\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='collect-results.py', description=__doc__)
	parser.add_argument('-files', metavar='FILES', nargs='+', required=True, help='All files to collect numbers from.')
	parser.add_argument('-metric', metavar='METRIC', required=True, choices=['concordance', 'vcfeval-typable', 'truvari-typable'], help='evaluation metric.')
	parser.add_argument('-outfile', metavar='OUTFILE', required=True, help='Name of the output file.')
	args = parser.parse_args()
	if args.metric == 'concordance':
		collect_concordances(args.files, args.outfile)
	if args.metric == 'vcfeval-typable':
		collect_precision_recall(args.files, args.outfile, "vcfeval")
	if args.metric == "truvari-typable":
		collect_precision_recall(args.files, args.outfile, "truvari")
