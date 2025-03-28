if config["consensus_haplotypes"]:
	# make sure that all fields are provided
	for callset in config["consensus_haplotypes"]:
		if not "agc" in config["consensus_haplotypes"][callset]:
			raise RuntimeError("No agc file provided for callset " + callset)
			sys.exit(1)
		if not "names" in config["consensus_haplotypes"][callset]:
			raise RuntimeError("No names file provided for callset " + callset)

if config["phased_vcfs"]:
	for callset in config["phased_vcfs"]:
		if not "vcf" in config["phased_vcfs"][callset]:
			raise RuntimeError("No phased vcf provided for callset " + callset)
		else:
			# make sure VCF is gzip compressed
			if not config["phased_vcfs"][callset]["vcf"].endswith(".vcf.gz"):
				raise RuntimeError("VCF must be gzip compressed for callset " + callset)
		if not "reference" in config["phased_vcfs"][callset]:
			raise RuntimeError("No reference provided for callset " + callset)
		else:
			# make sure ".fai" index exists
			if not os.path.isfile(config["phased_vcfs"][callset]["reference"] + ".fai"):
				raise RuntimeError("fai index missing for reference of callset " + callset)

if config["unphased_vcfs"]:
	for callset in config["unphased_vcfs"]:
		if not "vcf" in config["unphased_vcfs"][callset]:
			raise RuntimeError("No unphased vcf provided for callset " + callset)
		else:
			# make sure VCF is gzip compressed
			if not config["unphased_vcfs"][callset]["vcf"].endswith(".vcf.gz"):
				raise RuntimeError("VCF must be gzip compressed for callset " + callset) 
		if not "reference" in config["unphased_vcfs"][callset]:
			raise RuntimeError("No reference provided for callset " + callset)
		else:
			# make sure ".fai" index exists
			if not os.path.isfile(config["unphased_vcfs"][callset]["reference"] + ".fai"):
				raise RuntimeError("fai index missing for reference of callset " + callset)
		if not "sample_sheet" in config["unphased_vcfs"][callset]:
			raise RuntimeError("No sample sheet provided for callset " + callset)
		if not "genetic_maps" in config["unphased_vcfs"][callset]:
			raise RuntimeError("No genetic maps provided for callset " + callset)

