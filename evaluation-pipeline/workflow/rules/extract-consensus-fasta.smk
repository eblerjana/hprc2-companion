rule extract_from_agc:
	input:
		lambda wildcards: CONSENSUS_HAPLOTYPES[wildcards.callset]["agc"]
	output:
		consensus = "{results}/evaluation/{callset}/{sample}_{haplotype}/{sample}_{haplotype}_consensus.fa.gz",
		contig_names = "{results}/evaluation/{callset}/{sample}_{haplotype}/{sample}_{haplotype}_contigs.txt",
		temp_fa = temp("{results}/evaluation/{callset}/{sample}_{haplotype}/{sample}_{haplotype}_temp.fa.gz")
	conda:
		"../envs/agc.yaml"
	params:
		name = lambda wildcards: CONSENSUS_NAMES[wildcards.callset][wildcards.sample][wildcards.haplotype],
		chromosomes = "--exclude " + " ".join(SEX_CHROMOSOMES) if SEX_CHROMOSOMES else ""
	resources:
		mem_mb = 40000
	shell:
		"""
		agc getset {input} {params.name} | bgzip > {output.temp_fa}
		samtools faidx {output.temp_fa}
		zcat {output.temp_fa} | python3 workflow/scripts/get-contig-names.py {params.chromosomes} > {output.contig_names}
                samtools faidx -r {output.contig_names} {output.temp_fa} | bgzip > {output.consensus}
                samtools faidx {output.consensus}
		"""
