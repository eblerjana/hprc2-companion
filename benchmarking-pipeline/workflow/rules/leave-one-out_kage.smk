
rule kage_index:
	"""
	Run Kage indexing.
	"""
	input:
		reference = REFERENCE,
		vcf = "{results}/leave-one-out/input-panel/panel-{sample}_bi.vcf"
	output:
		temp("{results}/leave-one-out/kage/index-{sample}/index-{sample}.npz")
	log:
		"{results}/leave-one-out/kage/index-{sample}/index-{sample}.log"
	resources:
		mem_mb = 800000,
		walltime = "10:00:00"
	benchmark:
		"{results}/leave-one-out/kage/index-{sample}/index-{sample}.benchmark.txt"
	singularity:
		"workflow/container/kage.sif"
	threads:
		24
	shell:
		"""
		kage index -r {input.reference} -v {input.vcf} -o {output} -k 31 -t {threads}  &> {log}
		"""

rule kage_genotype:
	"""
	Run Kage genotyping.
	"""
	input:
		reads = lambda wildcards: ILLUMINA[wildcards.sample],
		index = "{results}/leave-one-out/kage/index-{sample}/index-{sample}.npz"
	output:
		"{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping_tmp.vcf"
	log:
		"{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping.log"
	resources:
		mem_mb = 50000,
		walltime = "02:00:00"
	benchmark:
		 "{results}/leave-one-out/kage/{sample}/{sample}_kage.benchmark.txt"	
	singularity:
		"workflow/container/kage.sif"
	params:
		index = "{results}/leave-one-out/kage/index-{sample}/index-{sample}"
	threads: 24
	shell:
		"""
		kage genotype -i {params.index} -r {input.reads} -t {threads} --average-coverage {AVG_ILLUMINA_COV} -k 31 -o {output} &> {log}
		"""

rule kage_postprocess:
	"""
	Kage does not keep track of the variant IDs in the INFO field.
	This rule adds them back in (assuming the variant representation is identical)
	"""
	input:
		kage = "{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping_tmp.vcf",
		vcf = "{results}/leave-one-out/input-panel/panel-{sample}_bi.vcf.gz"
	output:
		"{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping.vcf"
	shell:
		"""
		cat {input.kage} | python3 workflow/scripts/add-ids.py {input.vcf} > {output}
		"""
