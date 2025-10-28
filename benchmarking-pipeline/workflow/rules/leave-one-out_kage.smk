


rule kage_prepare_fasta:
	"""
	KAGE cannot handle extra strings after the chromosome name,
	even when separated by a whitespace. This rule removes such
	strings.
	"""
	input:
		REFERENCE
	output:
		"{results}/leave-one-out/kage/reference.fa"
	benchmark:
		 "{results}/leave-one-out/kage/reference.benchmark"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		cut -d" " -f1 {input} > {output}
		samtools faidx {output}
		"""



rule kage_extract_chromosome_from_vcf:
	"""
	Extract a chromosome from the VCF.
	"""
	input:
		"{results}/leave-one-out/input-panel/panel-{sample}_bi.vcf.gz"
	output:
		temp("{results}/leave-one-out/kage/panel-{sample}_{chrom}.vcf")
	benchmark:
		"{results}/leave-one-out/kage/panel-{sample}_{chrom}.benchmark"
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


rule kage_index:
	"""
	Run Kage indexing.
	"""
	input:
		reference = "{results}/leave-one-out/kage/reference.fa",
		vcf = "{results}/leave-one-out/kage/panel-{sample}_{chrom}.vcf"
	output:
		temp("{results}/leave-one-out/kage/index-{sample}/index-{sample}_{chrom}.npz")
	log:
		"{results}/leave-one-out/kage/index-{sample}/index-{sample}_{chrom}.log"
	resources:
		mem_mb = 500000,
		walltime = "05:00:00"
	wildcard_constraints:
		chrom = "|".join(KAGE_CHROMOSOMES)
	benchmark:
		"{results}/leave-one-out/kage/index-{sample}/index-{sample}_{chrom}.benchmark.txt"
	singularity:
		"workflow/container/kage.sif"
	threads:
		32
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
		index = "{results}/leave-one-out/kage/index-{sample}/index-{sample}_{chrom}.npz"
	output:
		temp("{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping_tmp_{chrom}.vcf")
	log:
		"{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping_{chrom}.log"
	resources:
		mem_mb = 50000,
		walltime = "05:00:00"
	wildcard_constraints:
		chrom = "|".join(KAGE_CHROMOSOMES)
	benchmark:
		 "{results}/leave-one-out/kage/{sample}/{sample}_kage_{chrom}.benchmark.txt"	
	singularity:
		"workflow/container/kage.sif"
	threads: 24
	shell:
		"""
		kage genotype -i {input.index} -r {input.reads} -t {threads} --average-coverage {AVG_ILLUMINA_COV} -k 31 -o {output} &> {log}
		"""

rule kage_postprocess:
	"""
	Kage does not keep track of the variant IDs in the INFO field.
	This rule adds them back in (assuming the variant representation is identical)
	"""
	input:
		kage = "{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping_tmp_{chrom}.vcf",
		vcf = "{results}/leave-one-out/input-panel/panel-{sample}_bi.vcf.gz"
	output:
		vcf = temp("{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping_{chrom}.vcf.gz"),
		tbi = temp("{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping_{chrom}.vcf.gz.tbi")
	resources:
		mem_mb = 50000
	shell:
		"""
		cat {input.kage} | python3 workflow/scripts/add-ids.py {input.vcf} | bgzip > {output}
		tabix -p vcf {output}
		"""


rule kage_concat_vcfs:
	"""
	Combine chromosome-specific VCFs into a single VCF.
	"""
	input:
		vcfs = expand("{{results}}/leave-one-out/kage/{{sample}}/{{sample}}_kage_bi_genotyping_{chrom}.vcf.gz", chrom = KAGE_CHROMOSOMES),
		tbis= expand("{{results}}/leave-one-out/kage/{{sample}}/{{sample}}_kage_bi_genotyping_{chrom}.vcf.gz.tbi", chrom = KAGE_CHROMOSOMES)
	output:
		"{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping.vcf"
	log:
		 "{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping.log"
	benchmark:
		 "{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping.benchmark"
	conda:
		"../envs/genotyping.yml"
	threads: 15
	resources:
		mem_mb = 50000,
		walltime = "05:00:00"
	shell:
		"""
		 bcftools concat -o {output} -Ov --threads {threads} {input.vcfs} &> {log}
		"""
