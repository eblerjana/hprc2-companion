import matplotlib.pyplot as plt
from collections import defaultdict
import sys
import math
import pandas as pd

plt.style.use('tableau-colorblind10')
plt.rcParams["font.family"] = "Nimbus Sans"


index = []
is_first = True
sample_to_changes = defaultdict(list)

for filename in sys.stdin:
	sample_name = filename.split('_')[-3]
	for line in open(filename.strip(), 'r'):
		fields = line.strip().split()
		if line.startswith("ref_HT"):
			continue
		sample_to_changes[sample_name].append(float(fields[1])) #  - float(fields[2]))
		if is_first:
			index.append(fields[0])
	is_first = False

print(sample_to_changes)

df = pd.DataFrame(sample_to_changes, index=index)
print(df)
df.plot.bar(rot=90)
plt.ylabel("change in variant-based QV")
plt.tight_layout()
plt.savefig(sys.argv[1])
