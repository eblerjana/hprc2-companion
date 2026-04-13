import gzip

CONSENSUS_HAPLOTYPES = config["consensus_haplotypes"]
PHASED_VCFS = config["phased_vcfs"]
UNPHASED_VCFS = config["unphased_vcfs"]
OUTNAME = config["outname"]
CONSENSUS_NAMES = {}
READS = {}

OVERLAPS = config["overlaps"] if "overlaps" in config else {}

# parse consensus names
for callset in CONSENSUS_HAPLOTYPES:
	CONSENSUS_NAMES[callset] = {}
	for line in open(CONSENSUS_HAPLOTYPES[callset]["names"], 'r'):
		fields = line.strip().split()
		if len(fields) < 3:
			raise RuntimeError("names file for callset " + callset + " is malformatted.")
			sys.exit(1)
		sample = fields[0]
		CONSENSUS_NAMES[callset][sample] = {}
		CONSENSUS_NAMES[callset][sample]["hap1"] = fields[1]
		CONSENSUS_NAMES[callset][sample]["hap2"] = fields[2]


# determine list of samples in VCFs
for callset in PHASED_VCFS:
	for line in gzip.open(PHASED_VCFS[callset]["vcf"], 'rt'):
		if line.startswith("##"):
			continue
		if line.startswith("#"):
			fields = line.strip().split()
			PHASED_VCFS[callset]["samples"] = fields[9:]
			break

for callset in UNPHASED_VCFS:
	for line in gzip.open(UNPHASED_VCFS[callset]["vcf"], 'rt'):
		if line.startswith("##"):
			continue
		if line.startswith("#"):
			fields = line.strip().split()
			UNPHASED_VCFS[callset]["samples"] = fields[9:]
			break


# add unphased callsets to PHASED_VCFS
for callset in UNPHASED_VCFS.keys():
	PHASED_VCFS[callset] = {}
	PHASED_VCFS[callset]["vcf"] = "{results}/phasing/{callset}_shapeit.bcf".format(callset = callset, results = OUTNAME)
	PHASED_VCFS[callset]["reference"] = UNPHASED_VCFS[callset]["reference"]
	PHASED_VCFS[callset]["samples"] = UNPHASED_VCFS[callset]["samples"]


# add newly generated consensus haplotypes
for callset in PHASED_VCFS:
	CONSENSUS_HAPLOTYPES[callset] = {}
	CONSENSUS_NAMES[callset] = {}
	CONSENSUS_HAPLOTYPES[callset]["agc"] = "{results}/haplotypes/{callset}.agc".format(callset = callset, results = OUTNAME)
	for sample in PHASED_VCFS[callset]["samples"]:
		CONSENSUS_NAMES[callset][sample] = {}
		CONSENSUS_NAMES[callset][sample]["hap1"] = callset + "_" + sample + "_hap1" 
		CONSENSUS_NAMES[callset][sample]["hap2"] = callset + "_" + sample + "_hap2"

# parse reads
for line in open(config['reads'], 'r'):
	fields = line.strip().split()
	if line.startswith('#sample'):
		continue
	sample = fields[0]
	READS[sample] = fields[1]

# for each callset, find overlap between reads and callset samples
SAMPLES = {}
for callset in CONSENSUS_NAMES:
	SAMPLES[callset] = []
	for sample in CONSENSUS_NAMES[callset].keys():
		if sample in READS.keys():
			SAMPLES[callset].append(sample)

SEX_CHROMOSOMES = config["sex_chromosomes"] if "sex_chromosomes" in config else []

#print(CONSENSUS_NAMES)
#print(CONSENSUS_HAPLOTYPES)
#print(PHASED_VCFS)
#print(UNPHASED_VCFS)
print(SAMPLES)
