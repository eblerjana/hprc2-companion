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

	metric = sys.argv[1]
	sample = sys.argv[2]
	outname = sys.argv[3]

	callset = None
	regions = []

	method_to_fscore = defaultdict(lambda: defaultdict(list))
	method_to_conc = defaultdict(lambda: defaultdict(list))
	method_to_precision_recall = defaultdict(lambda: defaultdict(list))

	method_to_name = {
		"kage-AoU": "KAGE-decomposed",
		"kage-AoU-bubble" : "KAGE-bubble",
		"pangenie-sampled-15-5x0.01" : "PanGenie v4"
	}

	region_to_name = {
		"autosomes" : "all",
		"easy-autosomes" : "easy regions",
		"repeat-autosomes": "repeat regions",
		"segdup-autosomes": "segdup regions"
	}


	for line in sys.stdin:
		if line.startswith("truthset"):
			continue
		fields = line.strip().split()
		m = fields[0]
		callset = fields[1]
		region = fields[2]
		filter = fields[3]
		method = fields[4]
		var = fields[5]
		s = fields[6]

		if not method in method_to_name:
			continue

		method = method_to_name[fields[4]]
		region = region_to_name[fields[2]]

		if m != metric:
			continue

		if s != sample:
			continue

		if region not in regions:
			regions.append(region)

		assert filter == "all"
		method_to_precision_recall[var][method + " (precision)"].append(float(fields[7]))
		method_to_precision_recall[var][method + " (recall)"].append(float(fields[8]))
		method_to_fscore[var][method].append(float(fields[9]))
		method_to_conc[var][method].append(float(fields[10]))




	with PdfPages(outname) as pdf:
		for variant in method_to_conc:
			bar_plot(dict(method_to_fscore[variant]), regions, "F-score (" + metric + ")", callset + ", " + sample + ", " + variant,  total_width=.8, single_width=.9)
			pdf.savefig()

			bar_plot(dict(method_to_precision_recall[variant]), regions, metric + " statistics", callset + ", " + sample + ", " + variant,  total_width=.8, single_width=.9, colors=["cornflowerblue", "royalblue", "sandybrown", "peru", "limegreen", "forestgreen"])
			pdf.savefig()


			if metric == "truvari":	
				bar_plot(dict(method_to_conc[variant]), regions, "Genotype concordance (truvari)", callset + ", " + sample + ", " + variant,  total_width=.8, single_width=.9)
				pdf.savefig()
