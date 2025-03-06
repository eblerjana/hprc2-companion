
#########################################################
##################     Evaluation     ###################
#########################################################


rule evaluation_prepare_truth:
	"""
	Extract ground truth genotypes for sample.
	"""
	input:
		"{results}/leave-one-out/pangenie/preprocessed-vcfs/{sample}_bi_no-missing.vcf.gz"
	output:
		temp("{results}/leave-one-out/pangenie/truth/{sample}-truth.vcf")
	conda:
		"../envs/genotyping.yml"
	priority: 1
	resources:
		mem_mb=20000
	log:
		"{results}/leave-one-out/pangenie/truth/{sample}-truth.log"
	shell:
		"bcftools view --samples {wildcards.sample} {input} 2> {log} 1> {output}"



rule evaluation_compress_vcf:
	input:
		"{results}/leave-one-out/{filename}.vcf"
	output:
		vcf = "{results}/leave-one-out/{filename}.vcf.gz",
		tbi = "{results}/leave-one-out/{filename}.vcf.gz.tbi"
	priority: 1
	shell:
		"""
		bgzip -c {input} > {output.vcf}
		tabix -p vcf {output.vcf}
		"""



rule evaluation_alleles_per_bubble:
	"""
	Compute allele statistics and regions of complex bubbles.
	"""
	input:
		PANEL_MULTI
	output:
		plot = "{results}/leave-one-out/alleles-per-bubble.pdf",
		bed = "{results}/leave-one-out/complex-bubbles.bed"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb=20000,
		walltime="1:00:00"
	shell:
		"""
		zcat {input} | python3 workflow/scripts/variant-statistics.py {output.plot} 1 > {output.bed}
		"""



rule evaluation_prepare_beds:
	"""
	Prepare a BED containing biallelic regions.
	"""
	input:
		bed = "{results}/leave-one-out/complex-bubbles.bed",
		fai = REFERENCE + '.fai'
	output:
		bed = "{results}/leave-one-out/biallelic-bubbles.bed",
		tmp = temp("{results}/leave-one-out/biallelic-bubbles.fai"),
		bed_tmp = temp("{results}/leave-one-out/biallelic-bubbles.tmp")
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		sort -k1,1d -k 2,2n -k 3,3n {input.fai} > {output.tmp}
		sort -k1,1d -k 2,2n -k 3,3n {input.bed} > {output.bed_tmp}
		bedtools complement -i {output.bed_tmp} -g {output.tmp} > {output.bed}
		"""



rule evaluation_untypable_ids:
	"""
	Determine untypable alleles based on unmerged callsets (i.e. independent of graph construction.
	This does not account for IDs that possible went missing during graph construction, i.e.
	the multi-allelic input VCF.
	"""
	input:
		PANEL_BI
	output:
		lists = "{results}/leave-one-out/pangenie/untypable-ids/{sample}-untypable.tsv",
		summary = temp("{results}/leave-one-out/pangenie/untypable-ids-{sample}.tsv")
	params:
		out = "{results}/leave-one-out/pangenie/untypable-ids"
	conda:
		"../envs/genotyping.yml"
	shell:
		"zcat {input} | python3 workflow/scripts/untypable-ids-single.py {params.out} {wildcards.sample} > {output.summary}"



rule evaluation_generally_untypable_ids:
	"""
	Determine IDs untypable because they have been filtered out during graph construction.
	Whenever a graph is constructed by merging callset into graph, some IDs might go missing (due to conflicts or ".").
	This step catches such cases as well. Since graph VCF might contain uncovered IDs, some untypable IDs might not get caught,
	that is why files are combined with other untypables (above).
	"""
	input:
		biallelic = PANEL_BI,
		multiallelic = "{results}/leave-one-out/pangenie/input-panel/panel-{sample}.vcf",
		untypables = "{results}/leave-one-out/pangenie/untypable-ids-{sample}.tsv"
	output:
		"{results}/leave-one-out/pangenie/untypable-ids/{sample}-untypable-all.tsv"
	resources:
		mem_mb = 50000
	log:
		"{results}/leave-one-out/pangenie/untypable-ids/{sample}-untypable-all.log"
	shell:
		"python3 workflow/scripts/untypable-ids-general.py {input.biallelic} {input.multiallelic} 2> {log} | cat - {input.untypables} | sort | uniq > {output}" 



rule evaluation_remove_untypable:
	"""
	Remove untypable variant alleles from VCF.
	"""
	input:
		vcf = "{results}/leave-one-out/pangenie/{path}{sample}{other}.vcf",
		ids = "{results}/leave-one-out/pangenie/untypable-ids/{sample}-untypable-all.tsv"
	output:
		vcf = "{results}/leave-one-out/pangenie/{path}{sample}{other}-typable-{vartype}.vcf.gz",
		tbi = "{results}/leave-one-out/pangenie/{path}{sample}{other}-typable-{vartype}.vcf.gz.tbi"
	wildcard_constraints:
		sample = "|".join(LEAVE_ONE_OUT_SAMPLES),
		vartype = "|".join(ALLOWED_VARIANTS)
	resources:
		mem_mb = 20000,
		walltime = "1:00:00"
	priority: 1
	shell:
		"""
		cat {input.vcf} | python3 workflow/scripts/skip-untypable.py {input.ids} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


rule evaluation_rtg_format:
	input:
		REFERENCE
	output:
		directory("{results}/leave-one-out/SDF")
	resources:
		mem_total_mb=20000
	conda:
		"../envs/genotyping.yml"
	priority: 1
	shell:
		"""
		rtg format -o {output} {input}
		"""



def region_to_bed(wildcards):
	if wildcards.regions == "biallelic":
		return "{results}/leave-one-out/biallelic-bubbles.bed".format(results=wildcards.results)
	if wildcards.regions == "multiallelic":
		return "{results}/leave-one-out/complex-bubbles.bed".format(results=wildcards.results)
	if wildcards.regions in EVALUATION_REGIONS:
		return EVALUATION_REGIONS[wildcards.regions]
	assert(False)



rule evaluation_vcfeval:
	"""
	Compute precision + recall statistics for evaluation.
	"""
	input:
		callset = "{results}/leave-one-out/pangenie/{mode}-{size}/{sample}/{sample}-pangenie_bi_genotyping-typable-{vartype}.vcf.gz",
		callset_tbi = "{results}/leave-one-out/pangenie/{mode}-{size}/{sample}/{sample}-pangenie_bi_genotyping-typable-{vartype}.vcf.gz.tbi",
		baseline = "{results}/leave-one-out/pangenie/truth/{sample}-truth-typable-{vartype}.vcf.gz",
		baseline_tbi = "{results}/leave-one-out/pangenie/truth/{sample}-truth-typable-{vartype}.vcf.gz.tbi",
		regions = region_to_bed,
		sdf = "{results}/leave-one-out/SDF"
	output:
		"{results}/leave-one-out/pangenie/{mode}-{size}/{sample}/precision-recall-typable/{regions}_{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	priority: 1
	wildcard_constraints:
		sample = "|".join(SAMPLES),
		vartype = "|".join(ALLOWED_VARIANTS)
	params:
		tmp = "{results}/leave-one-out/pangenie/{mode}-{size}/{sample}/precision-recall-typable/{regions}_{vartype}_tmp",
		outname = "{results}/leave-one-out/pangenie/{mode}-{size}/{sample}/precision-recall-typable/{regions}_{vartype}",
		which = "--all-records"
	resources:
		mem_mb = 30000,
		walltime = "1:40:00"
	shell:
		"""
		rtg vcfeval -b {input.baseline} -c {input.callset} -t {input.sdf} -o {params.tmp} --ref-overlap --evaluation-regions {input.regions} {params.which} --Xmax-length 30000 > {output}.tmp
		mv {params.tmp}/* {params.outname}/
		mv {output}.tmp {output}
		rm -r {params.tmp}
		"""



rule evaluation_collect_typed_variants:
	"""
	Determine variants that went into re-typing per category
	(needed for concordance computation).
	"""
	input:
		callset = "{results}/leave-one-out/pangenie/preprocessed-vcfs/{sample}_bi_no-missing.vcf.gz",
		regions= region_to_bed,
		ids="{results}/leave-one-out/pangenie/untypable-ids/{sample}-untypable-all.tsv"
	output:
		"{results}/leave-one-out/pangenie/genotyped-ids/{sample}_{regions}_{vartype}.tsv"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		sample = "|".join(LEAVE_ONE_OUT_SAMPLES),
		vartype = "|".join(ALLOWED_VARIANTS)
	resources:
		mem_mb=50000
	priority: 1
	shell:
		"zcat {input.callset} | python3 workflow/scripts/skip-untypable.py {input.ids} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bedtools intersect -header -a - -b {input.regions} -u -f 0.5 | python3 workflow/scripts/get_ids.py > {output}"



rule evaluation_genotype_concordances:
	"""
	Compute weighted genotype concordance for evaluation-
	"""
	input:
		callset = "{results}/leave-one-out/pangenie/{mode}-{size}/{sample}/{sample}-pangenie_bi_genotyping-typable-{vartype}.vcf.gz",
		callset_tbi = "{results}/leave-one-out/pangenie/{mode}-{size}/{sample}/{sample}-pangenie_bi_genotyping-typable-{vartype}.vcf.gz.tbi",
		baseline = "{results}/leave-one-out/pangenie/truth/{sample}-truth-typable-{vartype}.vcf.gz",
		baseline_tbi =  "{results}/leave-one-out/pangenie/truth/{sample}-truth-typable-{vartype}.vcf.gz.tbi",
		regions = region_to_bed,
		typed_ids = "{results}/leave-one-out/pangenie/genotyped-ids/{sample}_{regions}_{vartype}.tsv",
	output:
		tmp_vcf1 = temp("{results}/leave-one-out/pangenie/{mode}-{size}/{sample}/concordance/{regions}_{vartype}_base.vcf"),
		tmp_vcf2 = temp("{results}/leave-one-out/pangenie/{mode}-{size}/{sample}/concordance/{regions}_{vartype}_call.vcf"),
		summary = "{results}/leave-one-out/pangenie/{mode}-{size}/{sample}/concordance/{regions}_{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		sample = "|".join(LEAVE_ONE_OUT_SAMPLES),
		vartype = "|".join(ALLOWED_VARIANTS)
	log:
		"{results}/leave-one-out/pangenie/{mode}-{size}/{sample}/concordance/{regions}_{vartype}/summary.log"
	resources:
		mem_mb = 40000,
		walltime = "0:40:00"
	priority: 1
	shell:
		"""
		bedtools intersect -header -a {input.baseline} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf1}
		bedtools intersect -header -a {input.callset} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf2}
		python3 workflow/scripts/genotype-evaluation.py {output.tmp_vcf1} {output.tmp_vcf2} {input.typed_ids} --qual 0 2> {log} 1> {output.summary}
		"""


########################################################
##################     Plotting      ###################
########################################################



rule evaluation_collect_results:
	"""
	Collect statistics across all samples.
	"""
	input:
		lambda wildcards: expand("{{results}}/leave-one-out/pangenie/{{mode}}-{{size}}/{sample}/{{metric}}/{{regions}}_{{vartype}}/summary.txt", sample = LEAVE_ONE_OUT_SAMPLES),
	output:
		"{results}/leave-one-out/pangenie/plots/{metric}/{metric}_{mode}-{size}_{regions}_{vartype}.tsv"
	priority: 1
	shell:
		"""
		python3 workflow/scripts/collect-results.py -files {input} -metric {wildcards.metric} -outfile {output}
		"""


rule plot_across_samplings:
	"""
	Compare concordances across different sampling sizes for a fixed region.
	"""
	input:
		expand("{{results}}/leave-one-out/pangenie/plots/{{metric}}/{{metric}}_{mode}_{{regions}}_{vartype}.tsv", mode = ['sampled-' + s for s in SAMPLING_SIZES] + ['subset-' + s for s in SUBSETTING_SIZES], vartype = ALLOWED_VARIANTS)
	output:
		"{results}/leave-one-out/pangenie/plots/{metric}_{regions}.pdf"
	wildcard_constraints:
		metric = "concordance|precision-recall-typable"
	priority: 1
	conda:
		"../envs/genotyping.yml"
	params:
		sources = lambda wildcards: ' '.join(['sampled-' + str(s) + '_' + wildcards.regions for s in SAMPLING_SIZES] + ['subset-' + s for s in SUBSETTING_SIZES])
	shell:
		"""
		python3 workflow/scripts/plot-results.py -files {input} -outname {output} -metric {wildcards.metric} -sources {params.sources}
		"""


rule plot_concordance_vs_untyped:
	"""
	Compare concordance + untyped variants across different sampling sizes for a fixed region.
	"""
	input:
		expand("{{results}}/leave-one-out/pangenie/plots/concordance/concordance_{mode}_{{regions}}_{vartype}.tsv", mode = ['sampled-' + s for s in SAMPLING_SIZES] + ['subset-' + s for s in SUBSETTING_SIZES], vartype = ALLOWED_VARIANTS)
	output:
		"{results}/leave-one-out/pangenie/plots/concordance-vs-untyped_{regions}.pdf"
	priority: 1
	conda:
		"../envs/genotyping.yml"
	params:
		sources = lambda wildcards: ' '.join(['sampled-' + str(s) + '_' + wildcards.regions for s in SAMPLING_SIZES] + ['subset-' + s for s in SUBSETTING_SIZES])
	shell:
		"""
		python3 workflow/scripts/plot-results.py -files {input} -outname {output} -metric concordance-vs-untyped -sources {params.sources}
		"""


rule plot_resources:
	"""
	Plot resources (single core CPU time / max RSS)	for different sampling sizes.
	"""
	input:
		expand("{{results}}/leave-one-out/pangenie/{mode}/{sample}/{sample}-pangenie_multi_benchmark.txt", mode = ['sampled-' + s for s in SAMPLING_SIZES] + ['subset-' + s for s in SUBSETTING_SIZES], sample = LEAVE_ONE_OUT_SAMPLES)	
	output:
		"{results}/leave-one-out/pangenie/plots/resources.pdf"
	conda:
		"../envs/genotyping.yml"
	params:
		outname = "{results}/leave-one-out/pangenie/plots/resources"
	shell:
		"""
		python3 workflow/scripts/plot-resources.py -files {input} -outname {params.outname}
		"""
