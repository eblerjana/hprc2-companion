rule extract_from_agc:
	input:
		lambda wildcards: CONSENSUS_HAPLOTYPES[wildcards.callset][wildcards.region]["agc"]
	output:
		"{results}/evaluation/{callset}/{sample}_{region}_{haplotype}/{sample}_{region}_{haplotype}_consensus.fa.gz"
	conda:
		"../envs/agc.yaml"
	params:
		name = lambda wildcards: CONSENSUS_NAMES[wildcards.callset][wildcards.region][wildcards.sample][wildcards.haplotype]
	resources:
		mem_mb = 40000
	shell:
		"""
		agc getset {input} {params.name} | bgzip > {output}
		samtools faidx {output}
		"""
