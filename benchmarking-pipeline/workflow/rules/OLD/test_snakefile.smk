PANEL_MULTI = "test-data/MC-hgsvc3-hprc-chm13_filtered_ids.vcf.gz"
LEAVE_ONE_OUT_SAMPLES = ["NA18989"]
SAMPLING_SIZES = ["5", "10", "15"]

rule all:
	input:
		expand("{results}/leave-one-out/pangenie/plots/panel-stats.pdf", results = "results-test")
		

rule evaluation_panel_prepare_ground_truth:
	input:
		PANEL_MULTI	
	output:
		temp("{results}/leave-one-out/pangenie/truth/{sample}-truth_multi.vcf")
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/leave-one-out/pangenie/truth/{sample}-truth_multi.log"
	resources:
		mem_mb = 20000
	shell:
		"""
		bcftools view --samples {wildcards.sample} {input} 2> {log} 1> {output}
		"""


rule evaluation_panel_evaluate:
	"""
	Collect statistics on total number of ground truth alleles vs
	fraction of these alleles covered by sampled panel.
	"""
	input:
		panel = "test-data/{sample}-pangenie_multi_panel_{size}.vcf",
		truth = "{results}/leave-one-out/pangenie/truth/{sample}-truth_multi.vcf" 
	output:
		bubble_stats = "{results}/leave-one-out/pangenie/sampled-{size}/{sample}/panel-stats/panel_bubble_{sample}-{size}.tsv",
		summary_stats = "{results}/leave-one-out/pangenie/sampled-{size}/{sample}/panel-stats/panel_summary_{sample}-{size}.tsv"
	log:
		"{results}/leave-one-out/pangenie/sampled-{size}/{sample}/panel-stats/panel_{sample}-{size}.log"
	resources:
		mem_mb = 30000
	params:
		sampling_size = lambda wildcards: wildcards.size.split('-')[0]
	shell:
		"""
		python3 workflow/scripts/evaluate_panel.py -panel {input.panel} -truth {input.truth} -outname panel -size {params.sampling_size} -sample {wildcards.sample} &> {log}
		"""


rule evaluation_panel_plot:
	input:
		expand("{{results}}/leave-one-out/pangenie/sampled-{size}/{sample}/panel-stats/panel_summary_{sample}-{size}.tsv", sample = LEAVE_ONE_OUT_SAMPLES, size = SAMPLING_SIZES)
	output:
		tsv = temp("{results}/leave-one-out/pangenie/plots/panel-stats.tsv"),
		pdf = "{results}/leave-one-out/pangenie/plots/panel-stats.pdf"
	resources:
		mem_mb = 1000
	shell:
		"""
		awk 'FNR>1 || NR==1' {input} > {output.tsv}
		python3 workflow/scripts/plot_panel_stats.py -tsv {output.tsv} -outfile {output.pdf}
		"""
