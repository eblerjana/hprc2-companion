

####################################################################################################
#  find out which variants in the truth set are not contained in the pangenome graph (=untypables)
####################################################################################################


rule external_annotate_variants_callset:
	"""
	Assign a unique ID to each variant. Variant matching with panel will get the
	same ID assigned as in panel.
	"""
	input:
		callset = lambda wildcards: EXTERNAL[wildcards.truthset]["vcf"],
		panel = PANEL_BI
	output:
		"{results}/external-calls/evaluation/{truthset}/{evalsample}/truthset-{truthset}_{evalsample}_{vartype}.vcf.gz"
	resources:
		mem_total_mb=20000,
		runtime_hrs=1,
		runtime_min=59
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools norm -m -any {input.callset} | bcftools view --samples {wildcards.evalsample} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | python3 workflow/scripts/annotate.py {input.panel} | bgzip -c > {output}
		tabix -p vcf {output}
		"""


rule external_rtg_format_callsets:
	"""
	Create SDF file (required by truvari)
	"""
	input:
		reference = REFERENCE
	output:
		directory("{results}/external-calls/evaluation/SDF")
	resources:
		mem_total_mb = 20000
	conda:
		"../envs/genotyping.yml"
	shell:
		"rtg format -o {output} {input}"



rule external_determine_false_negatives_vcfeval:
	"""
	For each panel sample determine its false negatives comparing to the truth set.
	these are variants only in the truth, but not detected in the sample itself.
	Later, intersect all false negatives across the panel samples, to find out
	which variants are not present in the pangenome graph, but are in the truth.
	These variants are not accessible by re-genotyping methods, that are unable
	to detect variants themselves ("untypables").
	"""
	input:
		truthset = "{results}/external-calls/evaluation/{truthset}/{evalsample}/truthset-{truthset}_{evalsample}_snp-indel.vcf.gz",
		panel = PANEL_BI,
		reference = REFERENCE,
		ref_index = REFERENCE + '.fai',
		sdf = "{results}/external-calls/evaluation/SDF"
	output:
		sample_vcf = "{results}/external-calls/evaluation/{truthset}/{evalsample}/untypables-{evalsample}/samples/{sample}-vcfeval.vcf.gz",
		fn = "{results}/external-calls/evaluation/{truthset}/{evalsample}/untypables-{evalsample}/samples/{sample}/vcfeval/fn.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	params:
		tmp = "{results}/external-calls/evaluation/{truthset}/{evalsample}/untypables-{evalsample}/samples/{sample}/vcfeval_temp",
		outname = "{results}/external-calls/evaluation/{truthset}/{evalsample}/untypables-{evalsample}/samples/{sample}/vcfeval"
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=59
	log:
		"{results}/external-calls/evaluation/{truthset}/{evalsample}/untypables-{evalsample}/samples/{sample}/vcfeval.log"
	shell:
		"""
		bcftools view --samples {wildcards.sample} {input.panel} | bcftools view --min-ac 1 | python3 workflow/scripts/fix-header.py {input.ref_index} | bgzip -c > {output.sample_vcf}
		tabix -p vcf {output.sample_vcf}
		rtg vcfeval -b {input.truthset} -c {output.sample_vcf} -t {input.sdf} -o {params.tmp} --squash-ploidy  &> {log}
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""



rule external_determine_false_negatives_truvari:
	"""
	Same as previous rule, but for SVs using truvari instead of vcfeval.
	"""
	input:
		truthset = "{results}/external-calls/evaluation/{truthset}/{evalsample}/truthset-{truthset}_{evalsample}_sv.vcf.gz",
		panel = PANEL_BI,
		reference = REFERENCE,
		ref_index = REFERENCE + '.fai',
	output:
		sample_vcf = "{results}/external-calls/evaluation/{truthset}/{evalsample}/untypables-{evalsample}/samples/{sample}-truvari.vcf.gz",
		fn = "{results}/external-calls/evaluation/{truthset}/{evalsample}/untypables-{evalsample}/samples/{sample}/truvari/fn.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	params:
		tmp = "{results}/external-calls/evaluation/{truthset}/{evalsample}/untypables-{evalsample}/samples/truvari_temp",
		outname = "{results}/external-calls/evaluation/{truthset}/{evalsample}/untypables-{evalsample}/samples/{sample}/truvari"
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=59
	log:
		"{results}/external-calls/evaluation/{truthset}/{evalsample}/untypables-{evalsample}/samples/{sample}/truvari.log"
	shell:
		"""
		bcftools view --samples {wildcards.sample} {input.panel} | bcftools view --min-ac 1 | python3 workflow/scripts/fix-header.py {input.ref_index} | bgzip -c > {output.sample_vcf}
		tabix -p vcf {output.sample_vcf}
		truvari bench -b {input.truthset} -c {output.sample_vcf} -f {input.reference} -o {params.tmp} --pick ac -r 2000 -C 5000 --passonly --no-ref a  &> {log}
		bgzip {params.tmp}/fn.vcf
		tabix -p vcf {params.tmp}/fn.vcf.gz
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""


rule external_determine_unique:
	"""
	intersect the sets of FNs computed for each panel sample. The intersection then defines the set of unique/untypable variants
	"""
	input:
		expand("{{results}}/external-calls/evaluation/{{truthset}}/{{evalsample}}/untypables-{{evalsample}}/samples/{sample}/{{method}}/fn.vcf.gz", sample=PANEL_SAMPLES)
	output:	
		unique_tsv = "{results}/external-calls/evaluation/{truthset}/{evalsample}/untypables-{evalsample}/{truthset}-unique_{method}.tsv",
		unique_vcf = "{results}/external-calls/evaluation/{truthset}/{evalsample}/untypables-{evalsample}/{truthset}-unique_{method}.vcf"
	conda:
		"../envs/genotyping.yml"
	params:
		n_files = len(PANEL_SAMPLES)
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=30
	shell:
		"""
		bcftools isec -n={params.n_files} -w1 {input}  > {output.unique_vcf}
		grep -v '#' {output.unique_vcf} | cut -f 3  > {output.unique_tsv}
		"""


####################################################################################################
##  compare genotyping results to the truthset, excluding untypable variants which cannot be 
## genotyped correctly by a re-genotyper
####################################################################################################

rule external_extract_variant_type_callset:
	"""
	Remove untypables from VCF (callset or groundtruth) and extract only variants of a specific type
	"""
	input:
		vcf = lambda wildcards: "{results}/external-calls/evaluation/{truthset}/{evalsample}/truthset-{truthset}_{evalsample}_{vartype}.vcf.gz" if wildcards.set in EXTERNAL else "{results}/external-calls/{set}/{evalsample}/{evalsample}_{set}_bi_genotyping.vcf.gz",
		untypable = "{results}/external-calls/evaluation/{truthset}/{evalsample}/untypables-{evalsample}/{truthset}-unique_{method}.tsv",
		ref_index = REFERENCE + ".fai"
	output:
		vcf = "{results}/external-calls/evaluation/{truthset}/{evalsample}/{set}-typable_{evalsample}_{vartype}_{method}.vcf.gz",
		tbi = "{results}/external-calls/evaluation/{truthset}/{evalsample}/{set}-typable_{evalsample}_{vartype}_{method}.vcf.gz.tbi"
	wildcard_constraints:
		truthset = "|".join(EXTERNAL.keys()),
		set= "|".join([k for k in EXTERNAL.keys()] + GENOTYPERS),
		vartype="snp-indel|sv"
	resources:
		mem_total_mb=30000
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools view --samples {wildcards.evalsample} | python3 workflow/scripts/fix-header.py {input.ref_index} |  python3 workflow/scripts/skip-untypable.py {input.untypable} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


rule external_extract_variant_type_callset_all:
	"""
	Extract evaluation sample from VCF (callset or groundtruth) and extract only variants of a specific type
	"""
	input:
		vcf = lambda wildcards: "{results}/external-calls/evaluation/{truthset}/{evalsample}/truthset-{truthset}_{evalsample}_{vartype}.vcf.gz" if wildcards.set in EXTERNAL else "{results}/external-calls/{set}/{evalsample}/{evalsample}_{set}_bi_genotyping.vcf.gz",
		ref_index = REFERENCE + ".fai"
	output:
		vcf="{results}/external-calls/evaluation/{truthset}/{evalsample}/{set}-all_{evalsample}_{vartype}_{method}.vcf.gz",
		tbi="{results}/external-calls/evaluation/{truthset}/{evalsample}/{set}-all_{evalsample}_{vartype}_{method}.vcf.gz.tbi"
	wildcard_constraints:
		truthset = "|".join(EXTERNAL.keys()),
		set= "|".join([k for k in EXTERNAL.keys()] + GENOTYPERS),
		vartype="snp-indel|sv"
	resources:
		mem_total_mb=30000
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools view --samples {wildcards.evalsample} | python3 workflow/scripts/fix-header.py {input.ref_index} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


rule external_prepare_evaluation_beds:
	input:
		callable_regions = lambda wildcards: EXTERNAL[wildcards.truthset]["callable_regions"],
		bed = lambda wildcards: EVALUATION_REGIONS[wildcards.region] if wildcards.region != 'all' else EXTERNAL[wildcards.truthset]["callable_regions"]
	output:
		"{results}/external-calls/evaluation/{truthset}/bed-files/{truthset}_{region}.bed"
	resources:
		mem_total_mb=30000,
		runtime_hrs=1
	conda:
		"../envs/genotyping.yml"
	shell:
		"bedtools intersect -a {input.callable_regions} -b {input.bed} > {output}"


rule external_vcfeval_callsets:
	input:
		callset = "{results}/external-calls/evaluation/{truthset}/{evalsample}/{genotyper}-{filter}_{evalsample}_{vartype}_vcfeval.vcf.gz",
		callset_tbi = "{results}/external-calls/evaluation/{truthset}/{evalsample}/{genotyper}-{filter}_{evalsample}_{vartype}_vcfeval.vcf.gz.tbi",
		baseline = "{results}/external-calls/evaluation/{truthset}/{evalsample}/{truthset}-{filter}_{evalsample}_{vartype}_vcfeval.vcf.gz",
		baseline_tbi = "results/external-calls/evaluation/{truthset}/{evalsample}/{truthset}-{filter}_{evalsample}_{vartype}_vcfeval.vcf.gz.tbi",
		regions = "{results}/external-calls/evaluation/{truthset}/bed-files/{truthset}_{region}.bed",
		reference = REFERENCE,
		ref_index = REFERENCE + ".fai",
		sdf = "{results}/external-calls/evaluation/SDF"
	output:
		summary = "{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfeval-{evalsample}/vcfeval_{evalsample}_{genotyper}_{vartype}_{filter}_region-{region}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		vartype = "snp-indel",
		filter = "all|typable"
	log:
		"{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfeval-{evalsample}/vcfeval_{evalsample}_{genotyper}_{vartype}_{filter}_region-{region}.log"
	params:
		tmp = "{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfeval-{evalsample}/vcfeval_{evalsample}_{genotyper}_{vartype}_{filter}_region-{region}_tmp",
		outname = "{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfeval-{evalsample}/vcfeval_{evalsample}_{genotyper}_{vartype}_{filter}_region-{region}"
	resources:
		mem_total_mb=20000,
		runtime_hrs=0,
		runtime_min=40
	shell:
		"""
		rtg vcfeval -b {input.baseline} -c {input.callset} -t {input.sdf} -o {params.tmp} --evaluation-regions {input.regions} > {output.summary}.tmp
		mv {params.tmp}/* {params.outname}/
		mv {output.summary}.tmp {output.summary}
		rm -r {params.tmp}
		"""

rule external_truvari_callsets:
	input:
		callset = "{results}/external-calls/evaluation/{truthset}/{evalsample}/{genotyper}-{filter}_{evalsample}_{vartype}_vcfeval.vcf.gz",
		callset_tbi = "{results}/external-calls/evaluation/{truthset}/{evalsample}/{genotyper}-{filter}_{evalsample}_{vartype}_vcfeval.vcf.gz.tbi",
		baseline = "{results}/external-calls/evaluation/{truthset}/{evalsample}/{truthset}-{filter}_{evalsample}_{vartype}_vcfeval.vcf.gz",
		baseline_tbi = "{results}/external-calls/evaluation/{truthset}/{evalsample}/{truthset}-{filter}_{evalsample}_{vartype}_vcfeval.vcf.gz.tbi",
		regions = "{results}/external-calls/evaluation/{truthset}/bed-files/{truthset}_{region}.bed",
		reference = REFERENCE,
		ref_index = REFERENCE + ".fai"
	output:
		summary = "{results}/external-calls/evaluation/{truthset}/{evalsample}/truvari-{evalsample}/truvari_{evalsample}_{genotyper}_{vartype}_{filter}_region-{region}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		vartype = "sv",
		filter = "all|typable"
	log:
		"{results}/external-calls/evaluation/{truthset}/{evalsample}/truvari-{evalsample}/truvari_{evalsample}_{genotyper}_{vartype}_{filter}_region-{region}.log"
	params:
		tmp = "{results}/external-calls/evaluation/{truthset}/{evalsample}/truvari-{evalsample}/truvari_{evalsample}_{genotyper}_{vartype}_{filter}_region-{region}_tmp",
		outname = "{results}/external-calls/evaluation/{truthset}/{evalsample}/truvari-{evalsample}/truvari_{evalsample}_{genotyper}_{vartype}_{filter}_region-{region}",
		summary = "{results}/external-calls/evaluation/{truthset}/{evalsample}/truvari-{evalsample}/truvari_{evalsample}_{genotyper}_{vartype}_{filter}_region-{region}/summary.json"
	resources:
		mem_total_mb=20000,
		runtime_hrs=0,
		runtime_min=40
	shell:
		"""
		truvari bench -b {input.baseline} -c {input.callset} -f {input.reference} -o {params.tmp} --pick ac --passonly --includebed {input.regions} -r 2000 --no-ref a -C 5000 --refine &> {log}
		mv {params.tmp}/* {params.outname}/
		cp {params.summary} {output.summary}
		rm -r {params.tmp}
		"""


def input_external_collect_results(wildcards):
	filenames = []
	for vartype in EXTERNAL[wildcards.truthset]["variants"]:
		for sample in EXTERNAL[wildcards.truthset]["evaluation_samples"]:
			for genotyper in GENOTYPERS:
				for filter in ["all", "typable"]:
					if vartype == "snp-indel":
						filenames.append("{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfeval-{evalsample}/vcfeval_{evalsample}_{genotyper}_{vartype}_{filter}_region-{region}/summary.txt".format(results=wildcards.results, 
																													truthset = wildcards.truthset,
																													evalsample = sample,
																													genotyper = genotyper,
																													vartype = vartype,
																													filter = filter,
																													region = wildcards.region))

					else:
						assert vartype == "sv"
						filenames.append("{results}/external-calls/evaluation/{truthset}/{evalsample}/truvari-{evalsample}/truvari_{evalsample}_{genotyper}_{vartype}_{filter}_region-{region}/summary.txt".format(results=wildcards.results, 
                                                                                                                                                                                                                                        truthset = wildcards.truthset,
                                                                                                                                                                                                                                        evalsample = sample,
                                                                                                                                                                                                                                        genotyper = genotyper,
                                                                                                                                                                                                                                        vartype = vartype,
                                                                                                                                                                                                                                        filter = filter,
          	                                                                                                                                                                                                                        region = wildcards.region))
	return filenames


rule external_collect_results:
	input:
		input_external_collect_results
	output:
		tsv = "{results}/external-calls/evaluation/plots/{truthset}/external_{truthset}_{region}.tsv",
		pdf = "{results}/external-calls/evaluation/plots/{truthset}/external_{truthset}_{region}.pdf"
	conda:
		"../envs/genotyping.yml"
	params:
		outname = "{results}/external-calls/evaluation/plots/{truthset}/external"
	shell:
		"""
		ls {input} | python3 workflow/scripts/collect-external-stats.py -outname {params.outname} 
		"""
