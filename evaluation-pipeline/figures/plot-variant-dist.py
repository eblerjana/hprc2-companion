import sys
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

mapping = sys.argv[1]
fai = sys.argv[2]
outname = sys.argv[3]

plt.rcParams["font.family"] = "Nimbus Sans"

chrom_to_end = {}
for line in open(fai, 'r'):
	fields = line.strip().split()
	chrom_to_end[fields[0]] = int(fields[1])


region_to_color = {
	"AMPL": "skyblue",
	"PAR": "limegreen", 
	"XDR": "gold", 
	"XTR": "pink",
	"CEN": "maroon",
	"DYZ": "firebrick",
	"HET": "dimgray",
	"SAT": "blue",
	"TELO": "violet",
	"other": "lightgrey"
}

### background coloring of regions

chrom_to_intervals = defaultdict(list)
for bedfile in sys.stdin:
	for line in open(bedfile.strip(), 'r'):
		fields = line.strip().split()
		region = bedfile.strip().split(".")[-2].split('_')[-2]
		chrom_to_intervals[fields[0]].append((int(fields[1]), int(fields[2]), region))

### variant calls

chrom_to_positions = defaultdict(list)

for line in open(mapping, 'r'):
	if line.startswith("#"):
		continue
	fields = line.strip().split()
	chrom_to_positions[fields[0]].append(int(fields[1]))


with PdfPages(outname) as pdf:

	for chrom in chrom_to_intervals:
		split_points = []
		split_colors = []
		split_labels = []
		for r in sorted(chrom_to_intervals[chrom], key=lambda x: x[0]):
			split_points.append(r[0])
			split_points.append(r[1])
			split_colors.append(region_to_color["other"])
			split_colors.append(region_to_color[r[2]])
			split_labels.append("other")
			split_labels.append(r[2])		

		# add end coordinate of chromosome
		split_points.append(chrom_to_end[chrom])
		split_colors.append(region_to_color["other"])
		split_labels.append("other")


		print(chrom)
		print(split_points)
		print(split_colors)



		plt.figure()
		plt.title(chrom)

		## For each entry, draw the rectangle along with the color, transparency, zorder and label for the legend
		for i in range(len(split_points)):
			prev_point = 0 if i == 0 else split_points[i-1]
			plt.axvspan(prev_point, split_points[i], facecolor=split_colors[i], alpha=0.2, zorder=-10)

		patches = []
		for r,c in region_to_color.items():
			patch = mpatches.Patch(color=c, label=r)
			patches.append(patch)

		plt.legend(handles=patches)

		# plot histogram
		plt.hist(chrom_to_positions[chrom], bins=1000, color="black")
		plt.ylabel('Variant count')
		plt.xlabel('Position (bp)')

		pdf.savefig()
