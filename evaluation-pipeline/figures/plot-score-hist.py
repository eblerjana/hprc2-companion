import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import statistics
from collections import defaultdict

phred_scores = defaultdict(lambda: [])
nr_svs = defaultdict(lambda: [])

max_phred = 0
min_phred = float('inf')

max_nr = 0
min_nr = float('inf')

for filename in sys.stdin:
	for line in open(filename.strip(), 'r'):
		if line.startswith('ref_HT'):
			continue
		name = filename.strip().split('/')[-1].split('.')[0]
		fields = line.strip().split()
		if fields[16] != 'inf':
			phred = float(fields[16])
			svs = int(fields[18])
			phred_scores[name].append(phred)
			nr_svs[name].append(svs)

			max_phred = max(max_phred, phred)
			min_phred = min(min_phred, phred)
			max_nr = max(max_nr, svs)
			min_nr = min(min_nr, svs)

max_phred = int(max_phred)
min_phred = int(min_phred)
max_nr = int(max_nr)
min_nr = int(min_nr)

plt.style.use('tableau-colorblind10')

with PdfPages(sys.argv[1]) as pdf:

	# plot phred score historgram
	plt.figure()
	plt.title(sys.argv[2])

	for name in phred_scores:

		plt.hist(phred_scores[name], bins = range(min_phred - 1, max_phred + 2, 1), histtype='step')
		plt.axvline(statistics.median(phred_scores[name]), color='k', linestyle='dashed', linewidth=0.5)

		print('Phred scores for ' + name + ':' )
		print('min: ' + str(min(phred_scores[name])))
		print('max: ' + str(max(phred_scores[name])))
		print('mean: ' + str(statistics.mean(phred_scores[name])))
		print('median: ' + str(statistics.median(phred_scores[name])))

	plt.xlim([0,65])
#	plt.ylim([0,500])
	plt.xlabel('Variant-based QVs')
	plt.ylabel('Count')
	pdf.savefig()

	# plot histogram on the number of SVs by interval
	plt.figure()
	plt.title(sys.argv[2])

	for name in nr_svs:

		plt.hist(nr_svs[name], bins = range(min_nr - 1, max_nr + 2, 1), histtype='step')
		plt.axvline(statistics.median(nr_svs[name]), color='k', linestyle='dashed', linewidth=1)

		print('Number of variants per interval for ' + name + ':' )
		print('min: ' + str(min(nr_svs[name])))
		print('max: ' + str(max(nr_svs[name])))
		print('mean: ' + str(statistics.mean(nr_svs[name])))
		print('median: ' + str(statistics.median(nr_svs[name])))


	plt.xlim([0,150])
	plt.ylim([0,400])
	plt.xlabel('Number of variants (> 20bp)')
	plt.ylabel('Count')
	pdf.savefig()
