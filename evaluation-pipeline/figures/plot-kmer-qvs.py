import sys
import matplotlib.pyplot as plt
import numpy as np

samples = ["HG00609", "HG03009", "HG01258", "NA12877", "NA18983", "NA19331"]
sample_to_unpolished = {}
sample_to_polished = {}

haplogroup = {
	"HG00609": "O", 
	"HG03009": "H",
	"HG01258": "J",
	"NA12877": "R",
	"NA18983": "D",
	"NA19331": "E"
}

plt.style.use('tableau-colorblind10')
plt.rcParams["font.family"] = "Nimbus Sans"

for line in open(sys.argv[1], 'r'):
	fields = line.strip().split()
	sample_name = fields[0].split('_')[-3].split('-')[-1]
	qv = float(fields[3])
	sample_to_unpolished[sample_name] = qv

for line in open(sys.argv[2], 'r'):
	fields = line.strip().split()
	sample_name = fields[0].split('_')[-3].split('-')[-1]
	qv = float(fields[3])
	sample_to_polished[sample_name] = qv


fig, ax = plt.subplots(figsize=(4, 5))
x_vals = np.arange(len(samples))
ax.bar(x_vals, [sample_to_unpolished[s] for s in samples], 0.2, label="unpolished")
ax.bar(x_vals + 0.2, [sample_to_polished[s] for s in samples], 0.2, label="polished")
ax.set_xticks(x_vals + 0.1, [s + ' (' + haplogroup[s] + ')' for s in samples], rotation=90)
ax.set_title("Reconstruction of HG01596")
ax.set_ylabel("k-mer based QV")
ax.set_xlabel("polishing reference sequence")
ax.legend(loc="upper right", ncols=2)
plt.tight_layout()
plt.savefig(sys.argv[3])
