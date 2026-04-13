
rule merqury_extract_autosomes:
	"""
	Extract only autosomes.
	"""
	input:
		"{results}/evaluation/{callset}/{sample}_{haplotype}/{sample}_{haplotype}_consensus.fa.gz"
	output:
		consensus = temp("{results}/evaluation/qvs-kmers/temp/{callset}_{sample}_{haplotype}_consensus.fa.gz"),
		contig_names = "{results}/evaluation/qvs-kmers/{callset}_{sample}_{haplotype}_contigs.txt"
	conda:
		"../envs/shapeit.yaml"
	params:
		chromosomes = "--exclude " + " ".join(SEX_CHROMOSOMES) if SEX_CHROMOSOMES else ""
	shell:
		"""
		zcat {input} | python3 workflow/scripts/get-contig-names.py {params.chromosomes} > {output.contig_names}
		samtools faidx -r {output.contig_names} {input} | bgzip > {output.consensus}
		samtools faidx {output.consensus}
		"""
		

rule merqury_counting:
	"""
	Counting kmers using meryl
	"""
	input:
		lambda wildcards: READS[wildcards.sample]
	output:
		directory("{results}/evaluation/qvs-kmers/read-counts-meryl/{sample}.meryl")
	threads:
		24
	resources:
		mem_mb = 32768,
		mem_gb = 32,
		walltime = "05:00:00"
	conda:
		"../envs/merqury.yaml"
	log:
		"{results}/evaluation/qvs-kmers/read-counts-meryl/{sample}.log"
	shell:
		"""
		meryl count k=21 memory={resources.mem_gb} threads={threads} {input} output {output} &> {log}
		"""



rule merqury_evaluate_calls:
	input:
		counts = "{results}/evaluation/qvs-kmers/read-counts-meryl/{sample}.meryl", 
		assembly = "{results}/evaluation/qvs-kmers/temp/{callset}_{sample}_{haplotype}_consensus.fa.gz"
	output:
		qv = "{results}/evaluation/qvs-kmers/{callset}_{sample}_{haplotype}.qv",
		completeness = "{results}/evaluation/qvs-kmers/{callset}_{sample}_{haplotype}.completeness.stats"
	threads:
		24
	resources:
		mem_mb = 32768,
		walltime = "02:00:00",
		mem_gb = 32
	wildcard_constraints:
		haplotype = "hap1|hap2"
	conda:
		"../envs/merqury.yaml"
	log:
		"{results}/evaluation/qvs-kmers/{callset}_{sample}_{haplotype}.log"
	params:
		out_prefix = "{results}/evaluation/qvs-kmers/{callset}_{sample}_{haplotype}"
	shell:
		"""
		./workflow/scripts/compute-qv.sh {input.counts} {input.assembly} {params.out_prefix} {threads} {resources.mem_gb} &> {log}
		./workflow/scripts/compute-completeness.sh {input.counts} {input.assembly} {params.out_prefix} {threads} {resources.mem_gb} 
		"""


rule merqury_aggregate_qvs:
	input:
		lambda wildcards: expand("{{results}}/evaluation/qvs-kmers/{{callset}}_{sample}_{haplotype}.qv", sample = [s for s in SAMPLES[wildcards.callset] if s in READS], haplotype = ["hap1", "hap2"])
	output:
		"{results}/evaluation/qvs-kmers/{callset}-qvs-stats.tsv"
	shell:
		"""
		cat {input} > {output}
		"""


rule merqury_aggregate_completeness:
	input:
		lambda wildcards: expand("{{results}}/evaluation/qvs-kmers/{{callset}}_{sample}_{haplotype}.completeness.stats", sample = [s for s in SAMPLES[wildcards.callset] if s in READS], haplotype = ["hap1", "hap2"])
	output:
		"{results}/evaluation/qvs-kmers/{callset}-completeness-stats.tsv"
	shell:
		"""
		cat {input} > {output}
		"""

