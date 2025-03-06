
rule kage_index:
	input:
		reference = REFERENCE,
		vcf = "{results}/leave-one-out/input-panel/panel-{sample}_bi.vcf"
	output:
		"{results}/leave-one-out/kage/index-{sample}/index-{sample}.npz"
	log:
		"{results}/leave-one-out/kage/index-{sample}/index-{sample}.log"
	resources:
		mem_mb = 20000
	benchmark:
		"{results}/leave-one-out/kage/index-{sample}/index-{sample}.benchmark.txt"
	conda:
		"../envs/kage.yml"
	params:
		outname = "{results}/leave-one-out/kage/index-{sample}/index-{sample}"
	shell:
		"""
		kage index -r {input.reference} -v {input.vcf} -o {params.outname} -k 31 &> {log}
		"""

rule kage_genotype:
	input:
		reads = lambda wildcards: ILLUMINA[wildcards.sample],
		index = "{results}/leave-one-out/kage/index-{sample}/index-{sample}.npz"
	output:
		"{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping.vcf"
	log:
		"{results}/leave-one-out/kage/{sample}/{sample}_kage_bi_genotyping.log"
	resources:
		mem_mb = 10000
	benchmark:
		 "{results}/leave-one-out/kage/{sample}/{sample}_kage.benchmark.txt"	
	conda:
		"../envs/kage.yml"
	params:
		index = "{results}/leave-one-out/kage/index-{sample}/index-{sample}"
	threads: 24
	shell:
		"""
		kage genotype -i {params.index} -r {input.reads} -t {threads} --average-coverage {AVG_ILLUMINA_COV} -k 31 -o {output} &> {log}
		"""
