rule extract_from_agc:
	input:
		lambda wildcards: CONSENSUS_HAPLOTYPES[callset]["agc"]
	output:
		"{results}/evaluation/{callset}/{sample}_{haplotype}/{sample}_{haplotype}_consensus.fa.gz"
	conda:
		"../envs/agc.yaml"
	params:
		name = lambda wildcards: CONSENSUS_NAMES[wildcards.callset][wildcards.sample][wildcards.haplotype]
	resources:
		mem_mb = 20000
	shell:
		"""
		agc getset {input} {params.name} | bgzip > {output}
		samtools faidx {output}
		"""
