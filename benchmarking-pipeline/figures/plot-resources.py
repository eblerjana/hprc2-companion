import sys
from matplotlib.backends.backend_pdf import PdfPages
import gzip
import matplotlib.pyplot as plt
import matplotlib
import argparse
import numpy as np

def extract_resources(filename):
	cpu_time = 0.0
	max_rss = 0.0
	wallclock_time = 0.0
	for line in open(filename, 'r'):
		if line.startswith('s'):
			continue
		fields = line.strip().split()
		cpu_time = float(fields[-1])
		max_rss = float(fields[2]) * 0.001048576 # snakemake reports in MiB, convert to GB
		times = [f for f in fields[1].split(':')]
		wallclock_time += float(times[-1]) + 60.0 * float(times[-2])
		if len(times) > 2:
			wallclock_time += int(times[-3]) * 3600.0
	return wallclock_time, max_rss, cpu_time


def plot_resources(files, outname):
	plt.rcParams["font.family"] = "Nimbus Sans"
	colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#808080', '#f7ce37', '#bc202c', '#251188' ]
	runtimes = {}
	wallclock_times = {}
	rss = {}

	sample_to_name = {
		"pangenie-sampled-5-5x0.01": "PanGenie v4 5",
		"pangenie-sampled-10-5x0.01": "PanGenie v4 10",
		"pangenie-sampled-15-5x0.01": "PanGenie v4 15",
		"pangenie-sampled-25-5x0.01": "PanGenie v4 25",
		"pangenie-subset-108": "PanGenie v3"
	}

	samples = set([])
	sizes = set([])

	for file in files:
		sample = sample_to_name[file.strip().split('/')[-3]]
		size = file.strip().split('/')[-2]
		samples.add(sample)
		sizes.add(size)
		wallclock_time, max_rss, cpu_time = extract_resources(file)
		runtimes[(sample, size)] = cpu_time
		rss[(sample, size)] = max_rss
		wallclock_times[(sample, size)] = wallclock_time

	samples = ["PanGenie v4 5", "PanGenie v4 10", "PanGenie v4 15", "PanGenie v4 25", "PanGenie v3"]
	sizes = sorted(list(sizes))
	x_values = np.arange(len(samples))

	with PdfPages(outname +  '.pdf') as pdf:
		# plot runtimes
		width = 0.16
		fig, ax = plt.subplots(figsize=(3, 4.5))
		values = []
		for sample in samples:
			line_runtimes = [runtimes[(sample, size)] for size in sizes]
			values.append(line_runtimes)
		ax.boxplot(values, tick_labels=samples)
		ax.set_ylabel('Single core CPU seconds')
		ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
		plt.xticks(rotation=90)
		plt.tight_layout()
		pdf.savefig()
		plt.close()

		# plot wallclock times
		fig, ax = plt.subplots(figsize=(3, 4.5))
		values = []
		for sample in samples:
			line_wallclock = [wallclock_times[(sample, size)] for size in sizes]
			values.append(line_wallclock)
		ax.boxplot(values, tick_labels=samples)
		ax.set_ylabel('Wallclock seconds')
		ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
		plt.xticks(rotation=90)
		plt.tight_layout()
		pdf.savefig()
		plt.close()	

		# plot resources
		fig, ax = plt.subplots(figsize=(3, 4.5))
		values = []
		for sample in samples:
			line_rss = [ rss[(sample, size)] for size in sizes]
			values.append(line_rss)
		ax.boxplot(values, tick_labels=samples)
		ax.set_ylabel('Max RSS [GB]')
		ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
		plt.xticks(rotation=90)
		plt.tight_layout()
		pdf.savefig()
		plt.close()



parser = argparse.ArgumentParser(prog='plot-resources.py', description="Plot resources.")
parser.add_argument('-files', metavar='FILES', nargs='+', help='files with results per sample.')
parser.add_argument('-outname', metavar='OUTNAME', required=True, help='Name of the output file.')
args = parser.parse_args()

plot_resources(args.files, args.outname)
