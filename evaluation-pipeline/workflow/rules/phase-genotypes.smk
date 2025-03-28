

rule shapeit_extract_males:
	"""
	Collect names of all male samples.
	"""
	input:
		lambda wildcards: UNPHASED_VCFS[wildcards.callset]["sample_sheet"]
	output:
		"{results}/phasing/{callset}_haploid-samples.txt"
	shell:
		"awk '$5==1' {input} | cut -f 2 > {output}"



rule shapeit_extract_chromosome:
	"""
	Extract specific chromosome and set low quality genotypes
	to missing.
	"""
	input:
		lambda wildcards: UNPHASED_VCFS[wildcards.callset]["vcf"]
	output:
		"{results}/phasing/vcf/{callset}_{chrom}.vcf.gz"
	conda:
		"../envs/shapeit.yaml"
	log:
		"{results}/phasing/vcf/{callset}_{chrom}.log"
	benchmark:
		"{results}/phasing/vcf/{callset}_{chrom}.benchmark.txt"
	resources:
		mem_mb = 70000,
		walltime = "05:00:00"
	threads: 10
	shell:
		"""
		bcftools view -r {wildcards.chrom} {input} --threads {threads} | bcftools +setGT  -- -t q -n ./. -i 'FMT/GQ<10' 2> {log} | bgzip > {output}
		tabix -p vcf {output}
		"""

rule shapeit_prepare_trios:
	"""
	Prepare a file with trio information.
	"""
	input:
		lambda wildcards: UNPHASED_VCFS[wildcards.callset]["sample_sheet"]
	output:
		"{results}/phasing/{callset}_trios.ped"
	shell:
		"""
		awk '($3!=\"0\") && ($4!=\"0\")' {input} | cut -f2,3,4 > {output}
		"""


rule shapeit_phase_common:
	"""
	Phase a chromosome using shapeit.
	"""
	input:
		vcf = "{results}/phasing/vcf/{callset}_{chrom}.vcf.gz",
		fam = "{results}/phasing/{callset}_trios.ped",
		map = lambda wildcards: UNPHASED_VCFS[wildcards.callset]["genetic_maps"][wildcards.chrom],
		haploids = "{results}/phasing/{callset}_haploid-samples.txt"
	output:
		"{results}/phasing/{callset}_shapeit_{chrom}.bcf"
	log:
		"{results}/phasing/{callset}_shapeit_{chrom}.log"
	benchmark:
		"{results}/phasing/{callset}_shapeit_{chrom}-benchmark.txt"
	conda:
		"../envs/shapeit.yaml"
	threads: 32
        resources:
		mem_total_mb = 200000, # 100000
		runtime_hrs = 30
	params:
		haploids = lambda wildcards: "--haploids " + "{results}/phasing/{callset}_haploid-samples.txt".format(results = wildcards.results, callset = wildcards.callset)  if ("X" in wildcards.chrom) or ("Y" in wildcards.chrom) else ""
	shell:
		"""
		SHAPEIT5_phase_common --input {input.vcf} --pedigree {input.fam} --region {wildcards.chrom} {params.haploids} --map {input.map} --output {output} --thread {threads} &> {log}
		"""


rule shapeit_concat_vcfs:
	"""
	Combine the phased per-chromosome VCFs into a single one.
	"""
	input:
		lambda wildcards: expand("{{results}}/phasing/{{callset}}_shapeit_{chrom}.bcf", chrom = [c for c in UNPHASED_VCFS[wildcards.callset]["genetic_maps"].keys()])
	output:
		"{results}/phasing/{callset}_shapeit.bcf"
	conda:
		"../envs/shapeit.yaml"
	benchmark:
		"{results}/phasing/{callset}_shapeit.benchmark.txt"
	threads: 24
        resources:
		mem_mb = 100000,
		walltime = "05:00:00"
	log:
		"{results}/phasing/{callset}_shapeit.log"
	shell:
		"""
		bcftools concat -o {output} -O z --threads {threads} {input} &> {log}
		tabix -p vcf {output}
		"""

