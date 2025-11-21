import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse


def parse_precision_recall_vcfeval(filename):
	precision = None
	recall = None
	fscore = None
	for line in open(filename, 'r'):
		if 'None' in line:
			fields = line.strip().split()
			assert len(fields) == 8
			assert fields[0] == 'None'
			precision = float(fields[5])
			recall = float(fields[6])
			fscore = float(fields[7])
		if line.startswith('0 total baseline variants'):
			precision = 0.0
			recall = 0.0
			fscore = 0.0
	assert precision is not None
	assert recall is not None
	assert fscore is not None
	return precision, recall, fscore, 0.0


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
			precision = float(fields[-1][:-1])
		if "recall" in fields[0]:
			assert recall is None
			recall = float(fields[-1][:-1])
		if "f1" in fields[0]:
			assert fscore is None
			fscore = float(fields[-1][:-1])
		if "gt_concordance" in fields[0]:
			assert concordance is None
			concordance = float(fields[-1][:-1])
	assert precision is not None
	assert recall is not None
	assert fscore is not None
	return precision, recall, fscore, concordance


def print_data(data, truthset, region, filters, samples, vartypes, genotypers, outname):
	with open(outname, 'w') as outfile:
		outfile.write('\t'.join(['truthset', 'region', 'filter', 'genotyper', 'vartype', 'sample', 'precision', 'recall', 'fscore', 'gt_concordance']) + '\n')
		for filter in filters:
			for vartype in vartypes:
				for sample in samples:
					for genotyper in genotypers:
						precision = data[(truthset, region, filter, sample, vartype, genotyper)][0] 
						recall = data[(truthset, region, filter, sample, vartype, genotyper)][1]
						fscore = data[(truthset, region, filter, sample, vartype, genotyper)][2]
						concordance = data[(truthset, region, filter, sample, vartype, genotyper)][3]
						outfile.write('\t'.join([truthset, region, filter, genotyper, vartype, sample, str(precision), str(recall), str(fscore), str(concordance)]) + '\n')


	
def plot_data(data, truthset, region, filters, samples, vartypes, genotypers, outname):
	plt.style.use('tableau-colorblind10')
	plt.rcParams["font.family"] = "Nimbus Sans"
	var_to_name = {
		'snp-indel': 'SNPs + indels (1-49bp)',
		'sv': 'SVs (≥ 50bp)'
	}

	colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#808080', '#f7ce37', '#bc202c', '#251188']
	n_cols = len(vartypes)
	n_rows = len(filters) + 1 if 'all' in filters else len(filters)

	plt.figure(figsize=(20, 10))
	plot_index = 1
	for filter in filters:
		for var in vartypes:
			width = 0.16
			plt.subplot(n_rows, n_cols, plot_index)
			plot_index += 1
			for i,genotyper in enumerate(genotypers):
				fscores = []
				for sample in samples:
					fscores.append(data[(truthset, region, filter, sample, var, genotyper)][2])
				x_values = np.arange(len(samples))
				plt.bar(x_values + width*i, fscores, width=width, label=genotyper, color=colors[i])
			plt.title(var_to_name[var])
			plt.xticks(x_values + width*len(genotypers)/2 - width/2, samples)
			ylabel = ""
			if filter == "all":
				ylabel = "F-score"
			else:
				ylabel = "adjusted F-score"
			plt.ylabel(ylabel)
			plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.3)
			plt.tight_layout()
		# plot genotype concordances
		if filter == "all":
			for var in vartypes:
				width = 0.16
				plt.subplot(n_rows, n_cols, plot_index)
				plot_index += 1
				for i,genotyper in enumerate(genotypers):
					concordances = []
					for sample in samples:
						concordances.append(data[(truthset, region, filter, sample, var, genotyper)][3])
					x_values = np.arange(len(samples))
					plt.bar(x_values + width*i, concordances, width=width, label=genotyper, color=colors[i])
				plt.title(var_to_name[var])
				plt.xticks(x_values + width*len(genotypers)/2 - width/2, samples)
				ylabel = "gt_concordance"
				plt.ylabel(ylabel)
				plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.3)
				plt.tight_layout()
		

        # create legend
	handles = []
	labels = []
	for i, genotyper in enumerate(genotypers):
		label = genotyper
		line = matplotlib.lines.Line2D([],[], color=colors[i], markersize=100, linewidth=2.0, linestyle='-', label=label)
		handles.append(line)
		labels.append(label)
	plt.figlegend(handles, labels)
	plt.savefig(outname)



if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='collect-external-stats.py', description=__doc__)
#	parser.add_argument('-files', metavar='FILES', nargs='+', required=True, help='All files to collect numbers from.')
	parser.add_argument('-outname', metavar='OUTNAME', required=True, help='Prefix of output files.')
	args = parser.parse_args()

	truthsets = set([])
	regions = set([])
	filters = set([])
	samples = set([])
	vartypes = set([])
	genotypers = set([])

	data = {}

	for file in sys.stdin:
		truthset = file.strip().split('/')[-5]
		region = file.strip().split('/')[-2].split('_')[-1].split('region-')[-1]
		filter = file.strip().split('/')[-2].split('_')[-2]
		genotyper = file.strip().split('/')[-2].split('_')[-4]
		vartype = file.strip().split('/')[-2].split('_')[-3]
		sample = file.strip().split('/')[-4]
		if "vcfeval" in file:
			precision, recall, fscore, conc = parse_precision_recall_vcfeval(file.strip())
		else:
			precision, recall, fscore, conc = parse_precision_recall_truvari(file.strip())

		truthsets.add(truthset)
		regions.add(region)
		filters.add(filter)
		samples.add(sample)
		vartypes.add(vartype)
		genotypers.add(genotyper)

		data[(truthset, region, filter, sample, vartype, genotyper)] = [precision, recall, fscore, conc]

	truthsets = sorted(list(truthsets))
	regions = sorted(list(regions))
	filters = sorted(list(filters))
	samples = sorted(list(samples))
	vartypes = sorted(list(vartypes))
	genotypers = sorted(list(genotypers))


	for truthset in truthsets:
		for region in regions:
			print_data(data, truthset, region, filters, samples, vartypes, genotypers, args.outname + '_' + truthset + '_' + region + '.tsv')
			plot_data(data, truthset, region, filters, samples, vartypes, genotypers, args.outname + '_' + truthset + '_' + region + '.pdf')
