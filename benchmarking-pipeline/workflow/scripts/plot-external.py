import sys
from collections import defaultdict
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt


def bar_plot(data, xticks, ylabel, title, colors=None, total_width=0.8, single_width=1, legend=True):

	fig, ax = plt.subplots()

	if colors is None:
		colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

	# Number of bars per group
	n_bars = len(data)

	# The width of a single bar
	bar_width = total_width / n_bars

	# List containing handles for the drawn bars, used for the legend
	bars = []

	# Iterate over all data
	for i, (name, values) in enumerate(data.items()):
		# The offset in x direction of that bar
		x_offset = (i - n_bars / 2) * bar_width + bar_width / 2

		# Draw a bar for every value of that type
		for x, y in enumerate(values):
			bar = ax.bar(x + x_offset, y, width=bar_width * single_width, color=colors[i % len(colors)], label=[k for k in data.keys()][i])

	# Add a handle to the last drawn bar, which we'll need for the legend
	bars.append(bar[0])

	plt.title(title)
	ax.set_ylabel(ylabel)
	ax.set_xticks([i for i in range(len(values))], xticks)

	# Draw legend if we need
	if legend:

		handles = []
		labels = []
		for i, label in enumerate(data.keys()):
			line = matplotlib.lines.Line2D([],[], markersize=100, linewidth=2.0, linestyle='-', label=label, color=colors[i])
			handles.append(line)
			labels.append(label)
		plt.figlegend(handles, labels, loc='lower right', borderaxespad=5)


if __name__ == "__main__":

	plt.rcParams["font.family"] = "Nimbus Sans"

	vartype = sys.argv[1]
	sample = sys.argv[2]
	outname = sys.argv[3]

	callset = None

	regions = []

	method_to_fscore = defaultdict(list)
	method_to_ad_fscore = defaultdict(list)
	method_to_conc = defaultdict(list)

	for line in sys.stdin:
		if line.startswith("truthset"):
			continue
		fields = line.strip().split()
		callset = fields[0]
		region = fields[1]
		filter = fields[2]
		method = fields[3]
		var = fields[4]
		s = fields[5]

		if vartype != var:
			continue

		if s != sample:
			continue

		if region not in regions:
			regions.append(region)

		if filter == "all":
			method_to_fscore[method].append(float(fields[8]))
			method_to_conc[method].append(float(fields[9]))
		else:
			method_to_ad_fscore[method].append(float(fields[8]))


	with PdfPages(outname) as pdf:
		bar_plot(dict(method_to_fscore), regions, "F-score", callset + ", " + sample,  total_width=.8, single_width=.9)
		pdf.savefig()
		
		bar_plot(dict(method_to_ad_fscore), regions, "Adjusted F-score", callset + ", " + sample,  total_width=.8, single_width=.9)
		pdf.savefig()
	
		bar_plot(dict(method_to_conc), regions, "Genotype concordance (truvari)", callset + ", " + sample,  total_width=.8, single_width=.9)
		pdf.savefig()
