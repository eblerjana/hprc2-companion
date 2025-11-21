import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from collections import defaultdict

chrom_to_qvs = defaultdict(lambda: defaultdict(list))

for file in sys.stdin:
	for line in open(file.strip(), 'r'):
		if line.startswith('ref_HT'):
			continue
		fields = line.strip().split()
		chrom = fields[0].split('_')[-1].split(':')[0]
		if fields[16] == 'inf':
			chrom_to_qvs[chrom][file.strip()].append(-1)
		else:
			chrom_to_qvs[chrom][file.strip()].append(float(fields[16]))


	with PdfPages(sys.argv[1]) as pdf:
		for chrom in chrom_to_qvs:
			plt.figure()
			for file in chrom_to_qvs[chrom]:
				plt.title(chrom)
				y_values = chrom_to_qvs[chrom][file]
				x_values = [i*2 for i in range(len(y_values))]

				# ignore undefined QVs
				filtered_y = []
				filtered_x = []
				for x,y in zip(x_values, y_values):
					if y < 0:
						continue
					filtered_y.append(y)
					filtered_x.append(x)

				plt.bar(filtered_x, filtered_y, label = file.split('/')[-1])
			plt.ylabel('Window-wise QVs')
			plt.xticks([])
			plt.legend()
			pdf.savefig()
			plt.close()
