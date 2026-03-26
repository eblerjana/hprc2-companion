import sys
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages

name_to_outname = {
	"vienna_hap1": "ONT-S (hap1)",
	"vienna_hap2": "ONT-S (hap2)",
	"miller_hap1": "ONT-G (hap1)",
	"miller_hap2": "ONT-G (hap2)",
	"hifi_hap1": "HiFi (hap1)",
	"hifi_hap2": "HiFi (hap2)"
}

def get_name(name):
	if name in name_to_outname:
		return name_to_outname[name]
	else:
		return name

def get_variant_qv(fields):
	name = fields[0].split('_')[-2].split('-')[0] + " " +  get_name(fields[0].split('-')[-1])
	qv = float(fields[17])
	return name, qv

def get_kmer_qv(fields):
	name = fields[0].split('_')[-3].split('-')[0] + " " +  get_name(fields[0].split('-')[-1].split('_')[0] + '_' + fields[0].split('-')[-1].split('_')[1])
	qv = float(fields[3])
	return name, qv


qvs = defaultdict(dict)
names = []
first = True

outname = sys.argv[1]
qv_type = sys.argv[2]
assert qv_type in ['kmer', 'variant']
max_qv = 0

names = set([])

for filename in sys.stdin:

	if "unpolished" in filename:
		var_set = "unpolished"
	elif "onlySNP" in filename:
		var_set = "onlySNPs"
	elif "noSV" in filename:
		var_set = "noSVs"
	elif "repolish" in filename:
		var_set = "repolished"
	elif "all" in filename:
		var_set = "all"

	set_names = set([])
	for line in open(filename.strip(), 'r'):
		if line.startswith('ref'):
			continue
		fields = line.strip().split()
		name, qv = get_variant_qv(fields) if qv_type == 'variant' else get_kmer_qv(fields)
		qvs[var_set][name] = qv
		max_qv = max(qv, max_qv)
		set_names.add(name)

	if first:
		names = set_names
		first = False
	else:
		names = names.intersection(set_names)

names = sorted(list(names))

plt.style.use('tableau-colorblind10')
plt.rcParams["font.family"] = "Nimbus Sans"

with PdfPages(outname) as pdf:
	plt.figure()
	df = pd.DataFrame([	['unpolished'] + [qvs["unpolished"][n] for n in names],
			['SNPs'] + [qvs["onlySNPs"][n] for n in names],
			['SNPs,indels'] + [qvs["noSVs"][n] for n in names],
			['SNPs,indels,SVs'] + [qvs["all"][n] for n in names], 
			['repolished'] + [qvs["repolished"][n] for n in names]],
			columns = ['haplotypes'] + names  )
	ylabel = "Variant-based QV" if qv_type == 'variant' else "K-mer based QV"
	ax = df.plot(x='haplotypes',
		kind='bar',
		stacked=False, ylabel=ylabel, width=0.85, rot=0, xlabel="")

	for container in ax.containers:
		ax.bar_label(container,  fmt='%.2f', rotation=90)
	ax.set_ylim((0,max_qv + 7)) 

	plt.tight_layout()
	pdf.savefig()
	plt.close()

