configfile: "config/config.yaml"

CONSENSUS_HAPLOTYPES = config["consensus_haplotypes"]
PHASED_VCFS = config["phased_vcfs"]
UNPHASED_VCFS = config["unphased_vcfs"]

# determine list of samples in VCFs
for callset in PHASED_VCFS:
	for line in gzip.open(PHASED_VCFS[callset]["vcf"], 'rt'):
		if line.startswith("##"):
			continue
		if line.startswith("#"):
			PHASED_VCFS[callset]["samples"] = fields[9:]
			break

for callset in UNPHASED_VCFS:
	for line in gzip.open(UNPHASED_VCFS[callset]["vcf"], 'rt'):
		if line.startswith("##"):
			continue
		if line.startswith("#"):
			UNPHASED_VCFS[callset]["samples"] = fields[9:]
			break


# add unphased callsets to PHASED_VCFS
for callset in UNPHASED_VCFS.keys():
	PHASED_VCFS[callset] = {}
	PHASED_VCFS[callset]["vcf"] = "{results}/phasing/{callset}_shapeit.bcf".format(callset = callset)
	PHASED_VCFS[callset]["reference"] = UNPHASED_VCFS[callset]["reference"]
	PHASED_VCFS[callset]["samples"] = UNPHASED_VCFS[callset]["samples"]
