
rule genotyping_inputs:
	"""
	Helper rule.
	"""
	input:
		PANEL_MULTI,
		PANEL_BI
	output:
		"{results}/genotyping/input-panel-ready.txt"
	shell:
		"""
		echo "Panels ready" > {output}
		"""


rule genotyping_pangenie_prepare_panel:
	"""
	Prepare (uncompressed) PanGenie input panel.
	"""
	input:
		vcf = PANEL_MULTI
	output:
		temp("{results}/genotyping/panel-multi.vcf")
	shell:
		"""
		gunzip -c {input} > {output}
		"""


rule genotyping_pangenie_index:
	"""
	Create index for PanGenie.
	"""
	input:
		vcf = "{results}/genotyping/panel-multi.vcf",
		fasta = REFERENCE,
	output:
		temp(directory("{results}/genotyping/pangenie/index/"))
	log:
		"{results}/genotyping/pangenie/index.log"
	resources:
		mem_mb = 110000,
		walltime = "4:00:00"
	threads: 24
	params:
		out_prefix = "{results}/genotyping/pangenie/index/index"
	benchmark:
		"{results}/genotyping/pangenie/index.benchmark.txt"
	singularity:
		"workflow/container/pangenie.sif"
	shell:
		"""
		mkdir {output}
		PanGenie-index -v {input.vcf} -r {input.fasta} -o {params.out_prefix} -t {threads}  &> {log}
		"""


rule genotyping_pangenie_genotype_sampling:
	"""
	Run genotyping using sampling.
	"""
	input:
		reads = lambda wildcards: ILLUMINA[wildcards.sample],
		index = "{results}/genotyping/pangenie/index/"
	output:
		temp("{results}/genotyping/pangenie/{sample}_pangenie_multi_genotyping.vcf")
	log:
		"{results}/genotyping/pangenie/{sample}_pangenie_multi_genotyping.log"
	resources:
		mem_mb = 60000,
		walltime = "3:00:00"
	params:
		index = "{results}/genotyping/pangenie/index/index",
		out_prefix = "{results}/genotyping/pangenie/{sample}_pangenie_multi",
	benchmark:
		"{results}/genotyping/pangenie/{sample}/{sample}_pangenie.benchmark.txt"
	threads:
		24
	singularity:
		"workflow/container/pangenie.sif"
	shell:
		"""
		PanGenie -f {params.index} -i <(gunzip -c {input.reads}) -o {params.out_prefix} -t {threads} -j {threads} &> {log}
		"""


rule genotyping_convert_genotypes_to_biallelic:
	"""
	Convert genotyped VCF to biallelic representation.
	"""
	input:
		vcf = "{results}/genotyping/pangenie/{sample}_pangenie_multi_genotyping.vcf",
		biallelic = PANEL_BI
	output:
		bi = "{results}/genotyping/pangenie/{sample}_pangenie_bi_genotyping.vcf.gz",
		bi_tbi = "{results}/genotyping/pangenie/{sample}_pangenie_bi_genotyping.vcf.gz.tbi",
		multi = "{results}/genotyping/pangenie/{sample}_pangenie_multi_genotyping.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb=30000
	priority: 1
	shell:
		"""
		bgzip -c {input.vcf} > {output.multi}
		zcat {output.multi} | python3 workflow/scripts/convert-to-biallelic.py {input.biallelic} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' | bgzip > {output.bi}
		tabix -p vcf {output.bi}
		"""
