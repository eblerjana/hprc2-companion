

rule evaluation_panel_prepare_ground_truth:
	input:
		PANEL_MULTI	
	output:
		temp("{results}/leave-one-out/truth/{sample}-truth_multi.vcf")
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/leave-one-out/truth/{sample}-truth_multi.log"
	resources:
		mem_mb = 20000
	shell:
		"""
		bcftools view --samples {wildcards.sample} {input} 2> {log} 1> {output}
		"""


rule evaluation_panel_evaluate_stats:
	"""
	Collect statistics on total number of ground truth alleles vs
	fraction of these alleles covered by sampled panel.
	"""
	input:
		panel = "{results}/leave-one-out/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}_multi_panel.vcf.gz",
		truth = "{results}/leave-one-out/truth/{sample}-truth_multi.vcf" 
	output:
		bubble_stats = "{results}/leave-one-out/pangenie-sampled-{size}/{sample}/panel-stats/panel_bubble_{sample}_{size}.tsv",
		summary_stats = "{results}/leave-one-out/pangenie-sampled-{size}/{sample}/panel-stats/panel_summary_{sample}_{size}.tsv"
	log:
		"{results}/leave-one-out/pangenie-sampled-{size}/{sample}/panel-stats/panel_{sample}_{size}.log"
	resources:
		mem_mb = 30000
	conda:
		"../envs/pandas.yml"
	params:
		outname = "{results}/leave-one-out/pangenie-sampled-{size}/{sample}/panel-stats/panel"
	shell:
		"""
		python3 workflow/scripts/evaluate_panel.py -panel {input.panel} -truth {input.truth} -outname {params.outname} -size {wildcards.size} -sample {wildcards.sample} &> {log}
		"""


rule evaluation_panel_plot_missedAFs:
	"""
	Plot distribution of allele frequencies missing from the computed panel.
	"""
	input:
		expand("{{results}}/leave-one-out/pangenie-sampled-{size}/{sample}/panel-stats/panel_bubble_{sample}_{size}.tsv", sample = LEAVE_ONE_OUT_SAMPLES, size = SAMPLING_SIZES)
	output:
		"{results}/leave-one-out/plots/afs-missed-alleles.pdf"
	resources:
		mem_mb = 10000
	conda:
		"../envs/pandas.yml"
	shell:
		"""
		python3 workflow/scripts/plot_missed_afs.py -tsv {input} -outfile {output}
		"""


rule evaluation_panel_plot_stats:
	"""
	Plot panel statistics computed in previous rule across all samples + sampling_sizes.
	"""
	input:
		expand("{{results}}/leave-one-out/pangenie-sampled-{size}/{sample}/panel-stats/panel_summary_{sample}_{size}.tsv", sample = LEAVE_ONE_OUT_SAMPLES, size = SAMPLING_SIZES)
	output:
		tsv = temp("{results}/leave-one-out/plots/panel-stats.tsv"),
		pdf = "{results}/leave-one-out/plots/panel-stats.pdf"
	resources:
		mem_mb = 1000
	conda:
		"../envs/pandas.yml"
	shell:
		"""
		awk 'FNR>1 || NR==1' {input} > {output.tsv}
		python3 workflow/scripts/plot_panel_stats.py -tsv {output.tsv} -outfile {output.pdf}
		"""


rule evaluation_panel_collect_recomb:
	"""
	Collect chromsome-wise statistics on number of recombination events per haplotype
	of a specific sample.
	"""
	input:
		expand("{{results}}/leave-one-out/pangenie-sampled-{{size}}/{{sample}}/{{sample}}_pangenie-sampled-{{size}}_multi_paths_{chromosome}.tsv", chromosome = CHROMOSOMES)
	output:
		tsv = "{results}/leave-one-out/pangenie-sampled-{size}/{sample}/panel-stats/recombination_bubble_{sample}_{size}.tsv",
		summary_stats = "{results}/leave-one-out/pangenie-sampled-{size}/{sample}/panel-stats/recombination_summary_{sample}_{size}.tsv"
	params:
		outname = "{results}/leave-one-out/pangenie-sampled-{size}/{sample}/panel-stats/panel"
	resources:
		mem_mb = 30000
	conda:
		"../envs/pandas.yml"
	shell:
		"""
		awk 'FNR>1 || NR==1' {input} > {output.tsv}
		python3 workflow/scripts/evaluate_recombination.py -tsv {output.tsv} -size {wildcards.size} -sample {wildcards.sample} > {output.summary_stats}
		"""


rule evaluation_panel_plot_recombination:
	input:
		expand("{{results}}/leave-one-out/pangenie-sampled-{size}/{sample}/panel-stats/recombination_summary_{sample}_{size}.tsv", sample = LEAVE_ONE_OUT_SAMPLES, size = SAMPLING_SIZES)
	output:
		pdf = "{results}/leave-one-out/plots/recombination-stats.pdf"
	resources:
		mem_mb = 1000
	conda:
		"../envs/pandas.yml"
	params:
		chromosomes = " ".join([c for c in CHROMOSOMES])
	shell:
		"""
		python3 workflow/scripts/plot_recombination_stats.py -tsv {input} -outfile {output.pdf} -chromosomes {params.chromosomes}
		"""
