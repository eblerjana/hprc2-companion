
rule external_kage_prepare_panel:
	"""
	Prepare VCF by adding AC,AF,AN annotations.
	"""
	input:
		PANEL_BI
	output:
		vcf = temp("{results}/external-calls/panel-bi.vcf"),
		gz = temp("{results}/external-calls/panel-bi.vcf.gz")
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb = 20000
	shell:
		"""
		bcftools +fill-tags {input} -Oz -o {output.gz} -- -t AC,AN,AF
		tabix -p vcf {output.gz}
		gunzip -c {output.gz} > {output.vcf}
		"""


rule external_kage_prepare_fasta:
	"""
	KAGE cannot handle extra strings after the chromosome name,
	even when separated by a whitespace. This rule removes such
	strings.
	"""
	input:
		REFERENCE
	output:
		"{results}/external-calls/kage/reference.fa"
	benchmark:
		"{results}/external-calls/kage/reference.benchmark"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		cut -d" " -f1 {input} > {output}
		samtools faidx {output}
		"""


rule external_kage_extract_chromosome_from_vcf:
	"""
        Extract a chromosome from the VCF.
	"""
	input:
		"{results}/external-calls/panel-bi.vcf.gz"
	output:
		temp("{results}/external-calls/panel-bi_{chrom}.vcf")
	benchmark:
		"{results}/external-calls/panel-bi_{chrom}.benchmark"
	resources:
		mem_mb = 50000
	wildcard_constraints:
		chrom = "|".join(KAGE_CHROMOSOMES)
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools view -r {wildcards.chrom} {input} > {output}
		"""


rule external_kage_index:
	"""
	Run Kage indexing.
	"""
	input:
		reference = "{results}/external-calls/kage/reference.fa",
		vcf = "{results}/external-calls/panel-bi_{chrom}.vcf"
	output:
		temp("{results}/external-calls/kage/index/index_{chrom}.npz")
	log:
		"{results}/external-calls/kage/index/index_{chrom}.log"
	resources:
		mem_mb = 500000,
		walltime = "06:00:00"
	benchmark:
		"{results}/external-calls/kage/index/index_{chrom}.benchmark.txt"
	singularity:
		"workflow/container/kage.sif"
	wildcard_constraints:
		chrom = "|".join(KAGE_CHROMOSOMES)
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
		index = "{results}/external-calls/kage/index/index_{chrom}.npz"
	output:
		"{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping_tmp_{chrom}.vcf"
	log:
		"{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping_{chrom}.log"
	resources:
		mem_mb = 50000,
		walltime = "10:00:00"
	benchmark:
		 "{results}/external-calls/kage/{sample}/{sample}_kage_{chrom}.benchmark.txt"
	singularity:
		"workflow/container/kage.sif"
	wildcard_constraints:
		chrom = "|".join(KAGE_CHROMOSOMES)
	threads: 24
	shell:
		"""
		kage genotype -i {input.index} -r {input.reads} -t {threads} --average-coverage {AVG_ILLUMINA_COV} -k 31 -o {output} &> {log}
		"""


rule external_kage_postprocess:
	"""
	Kage does not keep track of the variant IDs in the INFO field.
	This rule adds them back in (assuming the variant representation is identical)
	"""
	input:
		kage = "{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping_tmp_{chrom}.vcf",
		vcf = "{results}/external-calls/panel-bi.vcf.gz"
	output:
		vcf = temp("{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping_{chrom}.vcf.gz"),
		tbi = temp("{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping_{chrom}.vcf.gz.tbi")
	resources:
		mem_mb = 50000,
		walltime = "05:00:00"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		cat {input.kage} | python3 workflow/scripts/add-ids.py {input.vcf} {wildcards.sample} | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


rule external_kage_concat_vcfs:
	"""
	Combine chromosome-specific VCFs into a single VCF.
	"""
	input:
		vcfs = expand("{{results}}/external-calls/kage/{{sample}}/{{sample}}_kage_bi_genotyping_{chrom}.vcf.gz", chrom = KAGE_CHROMOSOMES),
		tbis= expand("{{results}}/external-calls/kage/{{sample}}/{{sample}}_kage_bi_genotyping_{chrom}.vcf.gz.tbi", chrom = KAGE_CHROMOSOMES)
	output:
		"{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping.vcf.gz"
	log:
		"{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping.log"
	benchmark:
		"{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping.benchmark"
	conda:
		"../envs/genotyping.yml"
	threads: 15
	resources:
		mem_mb = 50000,
		walltime = "05:00:00"
	shell:
		"""
		bcftools concat -o {output} -Oz --threads {threads} {input.vcfs} &> {log}
		tabix -p vcf {output}
		"""

