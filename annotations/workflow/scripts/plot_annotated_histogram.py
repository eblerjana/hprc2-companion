import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from collections import defaultdict
import numpy as np

repeatnames = []

max_qv = 100
windowsize = 1000000
outname = sys.argv[1]

repeatclass_to_overlaps = defaultdict(lambda: [0] * max_qv)
qv_to_count = defaultdict(int)

for line in sys.stdin:
	fields = line.strip().split()
	if line.startswith('#'):
		repeatnames = fields[1:]
		continue

	if fields[3] == 'inf':
		continue

	qv = abs(int(float(fields[3])))
	qv_to_count[qv] += 1

	total_repeat_bp = 0

	for overlap, repeat in zip(fields[4:], repeatnames):
		bp_covered = float(overlap) * windowsize
		repeatclass_to_overlaps[repeat][qv] += bp_covered
		total_repeat_bp += bp_covered

	# compute bps in window not covered
	print(line)
	if total_repeat_bp > 1000000:
		assert total_repeat_bp < 1100000
		total_repeat_bp = 1000000
#	assert total_repeat_bp <= 1000000
	not_covered = 1000000 - total_repeat_bp

	repeatclass_to_overlaps["nonrepetitive"][qv] += not_covered

# normalize the counts by the windowsize
for repeatclass in repeatclass_to_overlaps:
	repeatclass_to_overlaps[repeatclass] = [r / windowsize for r in repeatclass_to_overlaps[repeatclass]]

# validate counts
computed_counts = [0] * max_qv
for l in repeatclass_to_overlaps.values():
	computed_counts = [sum(x) for x in zip(computed_counts, l)]
	
for qv in range(0, max_qv):
	print(qv_to_count[qv], computed_counts[qv])
	assert qv_to_count[qv] == round(computed_counts[qv])


# plot bars
x_values = [i for i in range(0, max_qv)]

plt.style.use('tableau-colorblind10')
plt.rcParams["font.family"] = "Nimbus Sans"
 
with PdfPages(outname) as pdf:
	fig, ax = plt.subplots()
	bottom = np.zeros(max_qv)
	for repeat, counts in repeatclass_to_overlaps.items():
		print(len(x_values))
		print(len(counts))
		print(len(bottom))
		p = ax.bar(x_values, counts, label=repeat, bottom=bottom)
		bottom += counts
	ax.legend(loc="upper right")
	ax.set_ylabel("Count")
	ax.set_xlabel("variant-based QV")
	pdf.savefig()
	plt.close()

