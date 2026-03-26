
####################################################################################################
##  Preprocessing
#####################################################################################################

rule external_normalize_truth:
	"""
	Normalize truth set and extract evaluation sample.
	"""
	input:
		callset = lambda wildcards: EXTERNAL[wildcards.truthset]["vcf"],
		ref_index = REFERENCE + '.fai',
		reference = REFERENCE
	output:
		"{results}/external-calls/evaluation/{truthset}/{evalsample}/truthset-{truthset}_{evalsample}_all.vcf.gz"
	resources:
		mem_mb=90000,
		walltime = "05:00:00"
	wildcard_constraints:
		truthset = "|".join(EXTERNAL.keys())
	conda:
		"../envs/bcftools.yml"
	shell:
		"""
		bcftools view --samples {wildcards.evalsample} {input.callset} | bcftools norm -m -any -f {input.reference} --check-ref x  | python3 workflow/scripts/fix-header.py {input.ref_index} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' | bgzip > {output}
		tabix -p vcf {output}
		"""


rule external_normalize_genotypes:
	"""
	Normalize the genotyped VCFs.
	"""
	input:
		vcf = "{results}/external-calls/{set}/{evalsample}/{evalsample}_{set}_bi_genotyping.vcf.gz",
		reference = REFERENCE,
		ref_index = REFERENCE + '.fai'
	output:
		"{results}/external-calls/evaluation/normalized-genotypes/{set}_{evalsample}_normalized_all.vcf.gz"
	resources:
		mem_mb = 50000,
		walltime = "02:00:00"
	conda:
		"../envs/bcftools.yml"
	shell:
		"""
		bcftools norm -f {input.reference} {input.vcf} | python3 workflow/scripts/fix-header.py {input.ref_index} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' | bgzip > {output}
		tabix -p vcf {output}
		"""


rule external_extract_variant_type:
	"""
	Extract only variants of a specific type (SNPs+indels / SVs)
	"""
	input:
		"{filename}_all.vcf.gz"
	output:
		"{filename}_{vartype}.vcf.gz"
	wildcard_constraints:
		vartype="snp-indel|sv"
	resources:
		mem_mb=30000
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		zcat {input} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output}
		tabix -p vcf {output}
		"""


rule external_rtg_format_callsets:
	"""
	Create SDF file (required by vcfeval)
	"""
	input:
		reference = REFERENCE
	output:
		directory("{results}/external-calls/evaluation/SDF")
	resources:
		mem_mb = 20000
	conda:
		"../envs/genotyping.yml"
	shell:
		"rtg format -o {output} {input}"


rule external_prepare_evaluation_beds:
	input:
		callable_regions = lambda wildcards: EXTERNAL[wildcards.truthset]["callable_regions"],
		bed = lambda wildcards: EVALUATION_REGIONS[wildcards.region] if wildcards.region != 'all' else EXTERNAL[wildcards.truthset]["callable_regions"]
	output:
		"{results}/external-calls/evaluation/{truthset}/bed-files/{truthset}_{region}.bed"
	resources:
		mem_mb=30000,
		walltime = "01:00:00"
	conda:
		"../envs/genotyping.yml"
	shell:
		"bedtools intersect -a {input.callable_regions} -b {input.bed} | bedtools sort > {output}"




####################################################################################################
##  Evaluation of small variants with vcfeval
####################################################################################################


rule external_vcfeval_callsets:
	input:
		callset = "{results}/external-calls/evaluation/normalized-genotypes/{set}_{evalsample}_normalized_snp-indel.vcf.gz",
		baseline = "{results}/external-calls/evaluation/{truthset}/{evalsample}/truthset-{truthset}_{evalsample}_snp-indel.vcf.gz",
		regions = "{results}/external-calls/evaluation/{truthset}/bed-files/{truthset}_{region}.bed",
		reference = REFERENCE,
		ref_index = REFERENCE + ".fai",
		sdf = "{results}/external-calls/evaluation/SDF"
	output:
		summary = "{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfeval-{evalsample}/vcfeval_{evalsample}_{set}_snp-indel_all_region-{region}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfeval-{evalsample}/vcfeval_{evalsample}_{set}_snp-indel_all_region-{region}.log"
	params:
		tmp = "{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfeval-{evalsample}/vcfeval_{evalsample}_{set}_snp-indel_all_region-{region}_tmp",
		outname = "{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfeval-{evalsample}/vcfeval_{evalsample}_{set}_snp-indel_all_region-{region}"
	resources:
		mem_mb=20000,
		walltime = "00:40:00"
	shell:
		"""
		rtg vcfeval -b {input.baseline} -c {input.callset} -t {input.sdf} -o {params.tmp} --evaluation-regions {input.regions} > {output.summary}.tmp
		mv {params.tmp}/* {params.outname}/
		mv {output.summary}.tmp {output.summary}
		rm -r {params.tmp}
		"""




####################################################################################################
###  Evaluation of SVs with truvari
####################################################################################################


rule external_truvari_callsets:
	input:
		callset = "{results}/external-calls/evaluation/normalized-genotypes/{set}_{evalsample}_normalized_sv.vcf.gz",
		baseline = "{results}/external-calls/evaluation/{truthset}/{evalsample}/truthset-{truthset}_{evalsample}_sv.vcf.gz",
		regions = "{results}/external-calls/evaluation/{truthset}/bed-files/{truthset}_{region}.bed",
		reference = REFERENCE,
		ref_index = REFERENCE + ".fai"
	output:
		summary = "{results}/external-calls/evaluation/{truthset}/{evalsample}/truvari-{evalsample}/truvari_{evalsample}_{set}_sv_all_region-{region}/summary.txt"
	conda:
		"../envs/truvari.yml"
	log:
		"{results}/external-calls/evaluation/{truthset}/{evalsample}/truvari-{evalsample}/truvari_{evalsample}_{set}_sv_all_region-{region}.log"
	params:
		tmp = "{results}/external-calls/evaluation/{truthset}/{evalsample}/truvari-{evalsample}/truvari_{evalsample}_{set}_sv_all_region-{region}_tmp",
		outname = "{results}/external-calls/evaluation/{truthset}/{evalsample}/truvari-{evalsample}/truvari_{evalsample}_{set}_sv_all_region-{region}",
		summary = "{results}/external-calls/evaluation/{truthset}/{evalsample}/truvari-{evalsample}/truvari_{evalsample}_{set}_sv_all_region-{region}/summary.json"
	resources:
		mem_mb=20000,
		walltime = "00:40:00"
	shell:
		"""
		truvari bench -b {input.baseline} -c {input.callset} -f {input.reference} -o {params.tmp} --pick multi -r 1000 -C 1000 -s 50 -S 15 --sizemax 100000 -p 0.0 -P 0.3 -O 0.0 --passonly --no-ref a --includebed {input.regions}  &> {log}
		mv {params.tmp}/* {params.outname}/
		cp {params.summary} {output.summary}
		rm -r {params.tmp}
		"""


####################################################################################################
###  Evaluation of all variants with vcfdist
####################################################################################################

rule external_naive_phasing:
	"""
	Fake phased genotypes by replacing / by |
	"""
	input:
		"{filename}.vcf.gz"
	output:
		"{filename}_phased.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		zcat {input} | sed '/^##/! s/\//\|/g' | bgzip > {output}
		tabix -p vcf {output}
		"""


rule external_vcfdist_callsets:
	input:
		callset = "{results}/external-calls/evaluation/normalized-genotypes/{set}_{evalsample}_normalized_all_phased.vcf.gz",
		baseline = "{results}/external-calls/evaluation/{truthset}/{evalsample}/truthset-{truthset}_{evalsample}_all_phased.vcf.gz",
		regions = "{results}/external-calls/evaluation/{truthset}/bed-files/{truthset}_{region}.bed",
		reference = REFERENCE,
		ref_index = REFERENCE + ".fai"
	output:
		summary = "{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfdist-{evalsample}/vcfdist_{evalsample}_{set}_all_all_region-{region}/precision-recall-summary.tsv"
	conda:
		"../envs/vcfdist.yml"
	log:
		"{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfdist-{evalsample}/vcfdist_{evalsample}_{set}_all_all_region-{region}.log"
	resources:
		mem_mb=600000,
		walltime = "10:00:00"
	threads: 24
	params:
		outname = "{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfdist-{evalsample}/vcfdist_{evalsample}_{set}_all_all_region-{region}/"
	shell:
		"""
		vcfdist {input.callset} {input.baseline} {input.reference}  --bed {input.regions} -p {params.outname} -t {threads} -r 200 &> {log}
		"""







def input_external_collect_results(wildcards):
	filenames = []
	for sample in EXTERNAL[wildcards.truthset]["evaluation_samples"]:
		for genotyper in GENOTYPERS:
			if wildcards.metric == "vcfeval":
				filenames.append("{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfeval-{evalsample}/vcfeval_{evalsample}_{genotyper}_snp-indel_all_region-{region}/summary.txt".format(
					results=wildcards.results,
					truthset = wildcards.truthset,
					evalsample = sample,
					genotyper = genotyper,
					region = wildcards.region))	
			elif wildcards.metric == "truvari":
				filenames.append("{results}/external-calls/evaluation/{truthset}/{evalsample}/truvari-{evalsample}/truvari_{evalsample}_{genotyper}_sv_all_region-{region}/summary.txt".format(
					results=wildcards.results,
					truthset = wildcards.truthset,
					evalsample = sample,
					genotyper = genotyper,
					region = wildcards.region))
			else:
				assert wildcards.metric == "vcfdist"	
				filenames.append("{results}/external-calls/evaluation/{truthset}/{evalsample}/vcfdist-{evalsample}/vcfdist_{evalsample}_{genotyper}_all_all_region-{region}/precision-recall-summary.tsv".format(
					results=wildcards.results,
					truthset = wildcards.truthset,
					evalsample = sample,
					genotyper = genotyper,
					region = wildcards.region))	
	return filenames


rule external_collect_results:
	input:
		input_external_collect_results
	output:
		tsv = "{results}/external-calls/evaluation/plots/{truthset}/external_{metric}_{truthset}_{region}.tsv",
		pdf = "{results}/external-calls/evaluation/plots/{truthset}/external_{metric}_{truthset}_{region}.pdf"
	conda:
		"../envs/genotyping.yml"
	params:
		outname = "{results}/external-calls/evaluation/plots/{truthset}/external_{metric}_{truthset}_{region}"
	shell:
		"""
		ls {input} | python3 workflow/scripts/collect-external-stats.py -outname {params.outname} 
		"""


rule external_plot_across_methods:
	input:
		expand("{{results}}/external-calls/evaluation/plots/{{truthset}}/external_{{metric}}_{{truthset}}_{region}.tsv", region = EVALUATION_REGIONS)
	output:
		"{results}/external-calls/evaluation/plots/{truthset}/{sample}_{metric}_{truthset}.pdf"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		cat {input} | python3 workflow/scripts/plot-external.py {wildcards.metric} {wildcards.sample} {output}
		"""
