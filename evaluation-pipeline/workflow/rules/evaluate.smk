rule evaluate_prepare_mapping:
	input:
		bed = lambda wildcards: ANNOTATION[wildcards.sample][wildcards.haplotype],
		vcf = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/{callset}_{sample}_{haplotype}.vcf.gz",
		fai = lambda wildcards: ASSEMBLIES[wildcards.sample][wildcards.haplotype] + ".fai"
	output:
		"{results}/evaluation/mapping_{callset}_{sample}_{haplotype}/{callset}_{sample}_{haplotype}_mapping.tsv"
	log:
		"{results}/evaluation/mapping_{callset}_{sample}_{haplotype}/{callset}_{sample}_{haplotype}_mapping.log"
	conda:
		"../envs/plotting.yaml"	
	params:
		outname = "{results}/evaluation/mapping_{callset}_{sample}_{haplotype}/{callset}_{sample}_{haplotype}"
	shell:
		"""
		python3 workflow/scripts/create-mapping.py {input.bed} {input.vcf} {input.fai} {params.outname} &> {log}
		"""
		


rule evaluate_calls:
	input:
		vcf = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/{callset}_{sample}_{haplotype}.vcf.gz",
		hap1 = "{results}/evaluation/{callset}/{sample}_hap1/{sample}_hap1_consensus.fa.gz",
		hap2 = "{results}/evaluation/{callset}/{sample}_hap2/{sample}_hap2_consensus.fa.gz",
		mapping = "{results}/evaluation/mapping_{callset}_{sample}_{haplotype}/{callset}_{sample}_{haplotype}_mapping.tsv"
	output:
		tsv = "{results}/evaluation/qvs/{callset}_{sample}_{haplotype}.tsv",
		bed_all = "{results}/evaluation/qvs/{callset}_{sample}_{haplotype}_all.bed",
		bed_err = "{results}/evaluation/qvs/{callset}_{sample}_{haplotype}_err.bed",
		intervals = "{results}/evaluation/qvs/{callset}_{sample}_{haplotype}_intervals.tsv"
	params:
		ref = "assembly"
	shell:
		"""
		zcat {input.vcf} | python3 workflow/scripts/evaluate-assemblies-intervals-mapping.py --name {wildcards.callset}_{wildcards.sample}_{wildcards.haplotype} --hap1 {input.hap1} --hap2 {input.hap2} --errors {output.bed_err} --all {output.bed_all} --reference {params.ref} --intervals {output.intervals} --mapping {input.mapping} > {output.tsv}
		"""


rule aggregate_results:
	input:
		lambda wildcards: expand("{{results}}/evaluation/qvs/{{callset}}_{sample}_{haplotype}.tsv", sample = SAMPLES[wildcards.callset], haplotype = ["hap1", "hap2"])
	output:
		"{results}/evaluation/qvs/{callset}-stats.tsv"
	shell:
		"""
		head -n 1 {input[0]} > {output}; tail -n +2 -q {input} >> {output}
		"""


rule plot_variants:
	input:
		"{results}/evaluation/qvs/{callset}_{sample}_{haplotype}_{which}.bed"
	output:
		"{results}/evaluation/qvs/plots/{callset}_{sample}_{haplotype}_{which}.pdf"
	wildcard_constraints:
		which = "all|err"
	conda:
		"../envs/plotting.yaml"
	params:
		name = "{callset}_{sample}_{haplotype}_{which}"
	resources:
		mem_mb = 10000	
	shell:
		"""
		cat {input} | python3 workflow/scripts/plot-svlen.py {params.name} {output}
		"""

rule plot_score_hist:
	input:
		"{results}/evaluation/qvs/{callset}_{sample}_{haplotype}_intervals.tsv"
	output:
		"{results}/evaluation/qvs/plots/{callset}_{sample}_{haplotype}_scores.pdf"
	conda:
		"../envs/plotting.yaml"
	log:
		"{results}/evaluation/qvs/plots/{callset}_{sample}_{haplotype}_scores.log"
	params:
		name = "{callset}_{sample}_{haplotype}"
	shell:
		"""
		cat {input} | python3 workflow/scripts/plot-score-hist.py {output} {params.name} &> {log}
		"""

#rule plot_scores_along_chromosome:
#	input:
#		"{results}/evaluation/qvs/{callset}_{sample}_{haplotype}_intervals.tsv"
#	output:
#		"{results}/evaluation/qvs/plots/{callset}_{sample}_{haplotype}_dist.pdf"
#	conda:
#		"../envs/plotting.yaml"
#	shell:
#		"""
#		"""
