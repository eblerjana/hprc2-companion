

rule filter_collect_samples:
	input:
		SAMPLE_SHEET
	output:
		unrelated = "{results}/filtering/unrelated-samples.tsv",
		all = "{results}/filtering/all-samples.tsv"
	run:
		with open(output.all, 'w') as all_samples, open(output.unrelated, 'w') as unrelated_samples, open(input[0], 'r') as infile:
			for line in infile:
				if line.startswith('#'):
					continue
				fields = line.strip().split()
				sample = fields[1]
				all_samples.write(sample + '\n')
				if (fields[2] == '0') and (fields[3] == '0'):
					# unrelated sample
					unrelated_samples.write(sample + '\n')
				else:
					# related sample
					related_samples.write(sample, '\n')


rule filter_extract_unrelated_samples:
	"""
	Extract only unrelated (parent) samples.
	"""
	input:
		vcf = "{results}/genotyping/pangenie_all-samples_unfiltered.vcf.gz",
		samples = "{results}/filtering/unrelated-samples.tsv"
	output:
		"{results}/genotyping/pangenie_unrelated-samples_unfiltered.vcf.gz"
	benchmark:
		"{results}/genotyping/pangenie_unrelated-samples_unfiltered.benchmark.txt"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb = 50000
	shell:
		"""
		bcftools view --sample-file {input.samples} --force-samples {input.vcf} -Oz {output}
		tabix -p vcf {output}
		"""


rule filter_extract_chromosome:
	"""
	Extract a chromosome from the VCF
	"""
	input:
		lambda wildcards: PANEL_BI if wildcards.callset == "panel" else "{results}/genotyping/pangenie_{subset}_unfiltered.vcf.gz"
	output:
		temp("{results}/filtering/chromosome-wise/{callset}_{chrom}.vcf.gz")
	benchmark:
		"{results}/filtering/chromosome-wise/{callset}_{chrom}.benchmark.txt"
	wildcard_constraints:
		callset = "panel|pangenie_all-samples_unfiltered|pangenie_unrelated-samples_unfiltered"
	resources:
		mem_mb = 50000
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools view -r {wildcards.chrom} {input} | bcftools +fill-tags -Oz -o {output} -- -t AN,AC,AF,AF_Hom,AF_Het
		"""


rule filter_compute_mendelian_consistency:
	"""
	Compute mendelian consistency.
	"""
	input:
		vcf = "{results}/filtering/chromosome-wise/pangenie_all-samples_unfiltered_{chrom}.vcf.gz",
		trios = SAMPLE_SHEET,
		samples = "{results}/filtering/all-samples.tsv"
	output:
		"{results}/filtering/pangenie_all-samples_unfiltered_mendelian-consistency_{chrom}.tsv"	
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/filtering/pangenie_all-samples_unfiltered_mendelian-consistency_{chrom}.log"
	benchmark:
		"{results}/filtering/pangenie_all-samples_unfiltered_mendelian-consistency_{chrom}.benchmark.tsv"
	resources:
		mem_mb = 30000,
		walltime = "08:00:00"
	shell:
		"""
		python3 workflow/scripts/evaluate-mendelian-consistency.py -vcf {input.vcf} -ped {input.trios} -samples {input.samples} -table {output} -column-prefix pangenie &> {log}
		"""


rule filter_compute_genotype_statistics:
	"""
	Compute genotype statistics.
	"""
	input:
		"{results}/filtering/chromosome-wise/{callset}_{chrom}.vcf.gz"
	output:
		"{results}/filtering/{callset}_genotype-statistics_{chrom}.tsv"
	wildcard_constraints:
		callset = "panel|pangenie_all-samples_unfiltered|pangenie_unrelated-samples_unfiltered"
	benchmark:
		"{results}/filtering/{callset}_genotype-statistics_{chrom}.benchmark.txt"
	resources:
		mem_mb = 30000,
		walltime = "08:00:00"
	params:
		flag = lambda wildcards: "" if wildcards.callset == "panel" else "--genotyping-stats"
	shell:
		"""
		python3 workflow/scripts/collect-vcf-stats.py {params.flag}
		"""	
		
			
