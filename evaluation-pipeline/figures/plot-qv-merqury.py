import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import argparse
from collections import defaultdict
import statistics

parser = argparse.ArgumentParser(prog='plot-qv.py', description=__doc__)
parser.add_argument('outname', metavar='OUTNAME', help='Name of the output PDF.')
parser.add_argument('-assembly', metavar='ASSEMBLYQVs', required=False, default='')
args = parser.parse_args()


plt.style.use('tableau-colorblind10')
plt.rcParams["font.family"] = "Nimbus Sans"

samples = set([])
qv_values = defaultdict(lambda: {})
callset_samples = set([])
assembly_samples = set([])

skip_samples = ["HG01890", "HG00731"]

for filename in sys.stdin:
	for line in open(filename.strip(), 'r'):
		fields = line.strip().split()
		file_info = fields[0].split('/')[-1].split('_')
		callset = file_info[0]
		sample = file_info[1]
		haplotype = file_info[2]

		assert haplotype in ['hap1', 'hap2']

		if sample in skip_samples:
			continue

		samples.add(sample)
		callset_samples.add(sample)
		qv_values[callset][sample + '_' + haplotype] = float(fields[3])



if args.assembly:
	qv_values["HGSVC3-assemblies"] = {}
	for line in open(args.assembly, 'r'):
		fields = line.strip().split()
		if line.startswith('hap'):
			continue
		if fields[0] in ['H0', 'wg']:
			continue
		sample = fields[-2].split('.')[0]
		samples.add(sample)
		assembly_samples.add(sample)
		assert fields[0] in ['H1', 'H2']
		haplotype = 'hap1' if fields[0] == 'H1' else 'hap2'
		qv_values["HGSVC3-assemblies"][sample + '_' + haplotype] = float(fields[-1])



y_values = defaultdict(lambda: [])
labels = []

# make sure all callset samples overlap with assembly samples
for s in callset_samples:
	assert s in assembly_samples

ordered_samples = sorted(callset_samples)
print(ordered_samples)

for sample in ordered_samples:
	labels.append(sample + '_hap1')
	labels.append(sample + '_hap2')
	for callset in qv_values.keys():
		for haplotype in ['hap1', 'hap2']:
			key_str = sample + '_' + haplotype
			value = qv_values[callset][key_str] if key_str in qv_values[callset] else None
			y_values[callset].append(value)

	
with PdfPages(args.outname) as pdf:
	x = [i*2 for i in range(len(labels))]

	plt.figure(figsize=(6,5))
	values = [ [y for y in y_values[k] if y is not None] for k in sorted(y_values.keys()) ]
	x_labels = [k for k in sorted(y_values.keys())]

	plt.boxplot(values)
	plt.xticks([i+1 for i in range(len(x_labels))], x_labels)
	plt.ylabel('k-mer based QV')
	plt.tight_layout()
	pdf.savefig()
	plt.close()


	plt.figure(figsize=(30,10))
	for callset in y_values.keys():
		plt.plot(x, y_values[callset], marker = 'o', label = callset)
	plt.xticks(x, labels, rotation='vertical')
	plt.legend()
	plt.ylabel('k-mer based QV')
	plt.tight_layout()
	pdf.savefig()
	plt.close()

	for v,k in zip(values, x_labels):
		print(k)
		print("Median: " +  str(statistics.median(v)))
		print("Min: " + str(min(v)))
		print("Max: " + str(max(v)))






