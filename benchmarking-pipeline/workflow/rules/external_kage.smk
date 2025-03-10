
rule external_kage_prepare_panel:
	input:
		vcf = PANEL_BI
	output:
		temp("{results}/external-calls/panel-bi.vcf")
	shell:
		"""
		gunzip -c {input} > {output}
		"""


rule external_kage_index:
	"""
	Run Kage indexing.
	"""
	input:
		reference = REFERENCE,
		vcf = "{results}/external-calls/panel-bi.vcf"
	output:
		temp("{results}/external-calls/kage/index/index.npz")
	log:
		"{results}/external-calls/kage/index/index.log"
	resources:
		mem_mb = 800000,
		walltime = "10:00:00"
	benchmark:
		"{results}/external-calls/kage/index/index.benchmark.txt"
	singularity:
		"workflow/container/kage.sif"
	threads:
		24
	shell:
		"""
		kage index -r {input.reference} -v {input.vcf} -o {output} -k 31 -t {threads}  &> {log}
		"""

rule external_kage_genotype:
	"""
	Run Kage genotyping.
	"""
	input:
		reads = lambda wildcards: ILLUMINA[wildcards.sample],
		index = "{results}/external-calls/kage/index/index.npz"
	output:
		"{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping_tmp.vcf"
	log:
		"{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping.log"
	resources:
		mem_mb = 50000,
		walltime = "02:00:00"
	benchmark:
		 "{results}/external-calls/kage/{sample}/{sample}_kage.benchmark.txt"	
	singularity:
		"workflow/container/kage.sif"
	params:
		index = "{results}/external-calls/kage/index/index"
	threads: 24
	shell:
		"""
		kage genotype -i {params.index} -r {input.reads} -t {threads} --average-coverage {AVG_ILLUMINA_COV} -k 31 -o {output} &> {log}
		"""

rule external_kage_postprocess:
	"""
	Kage does not keep track of the variant IDs in the INFO field.
	This rule adds them back in (assuming the variant representation is identical)
	"""
	input:
		kage = "{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping_tmp.vcf",
		vcf = "{results}/external-calls/input-panel/panel-{sample}_bi.vcf.gz"
	output:
		"{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping.vcf"
	shell:
		"""
		cat {input.kage} | python3 workflow/scripts/add-ids.py {input.vcf} > {output}
		"""
