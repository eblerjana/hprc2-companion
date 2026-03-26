import matplotlib.pyplot as plt
from collections import defaultdict
import sys
import math
import pandas as pd

plt.style.use('tableau-colorblind10')
plt.rcParams["font.family"] = "Nimbus Sans"

callsets = []

class_to_qvs = defaultdict(lambda: [])

for filename in sys.stdin:
	callset = filename.split('/')[-1].split('_')[0]
	callsets.append(callset)
	for line in open(filename.strip(), 'r'):
		if line.startswith("ref"):
			continue
		fields = line.strip().split()
		# combine all values computed across different chromosomes for the same class into one category.
		interval_id = fields[0]
		qv = float(fields[16])
		if math.isinf(qv):
			continue
		class_to_qvs[interval_id].append(qv)

print(callsets)

values = []
for c,v in class_to_qvs.items():
	values.append([c] + [qv for qv in v])

df = pd.DataFrame(values, columns = ['Region'] + callsets)
df.plot(x='Region',
        kind='bar',
        stacked=False,
	ylim=(-0.5, 55),
        title='')

plt.ylabel("variant-based QV")
plt.tight_layout()
plt.savefig(sys.argv[1])
