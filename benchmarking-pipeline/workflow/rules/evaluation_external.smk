

####################################################################################################
#  find out which variants in the truth set are not contained in the pangenome graph (=untypables)
####################################################################################################


rule annotate_variants_callset:
	"""
	Assign a unique ID to each variant. Variant matching with panel will get the
	same ID assigned as in panel.
	"""
	input:
		callset = lambda wildcards: EXTERNAL[wildcards.truthset]["vcf"]
		panel = PANEL_BI
		ref_index = REFERENCE + '.fai'
	output:
		"results/external-calls/evaluation/{truthset}/truthset-{truthset}_{vartype}.vcf.gz"
	resources:
		mem_total_mb=20000,
		runtime_hrs=1,
		runtime_min=59
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools norm -m -any {input.callset} | python3 workflow/scripts/fix-header.py {input.ref_index} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | python3 workflow/scripts/annotate.py {input.panel} | bgzip -c > {output}
		tabix -p vcf {output}
		"""


rule rtg_format_callsets:
	"""
	Create SDF file (required by truvari)
	"""
	input:
		reference = REFERENCE
	output:
		directory("results/external-calls/evaluation/SDF")
	resources:
		mem_total_mb = 20000
	conda:
		"../envs/genotyping.yml"
	shell:
		"rtg format -o {output} {input}"




# for each panel sample determine its false negatives comparing to the truth set.
# these are variants only in the truth, but not detected in the sample itself.
# later, intersect all false negatives across the panel samples, to find out
# which variants are not present in the pangenome graph, but are in the truth.
# these variants are not accessible by re-genotyping methods, that are unable
# to detect variants themselves ("untypables").
rule determine_false_negatives_vcfeval:
	input:
		truthset = "results/external-calls/evaluation/{truthset}/truthset-{truthset}_snp-indel.vcf.gz",
		panel = PANEL_BI,
		reference = REFERENCE,
		ref_index = REFERENCE + '.fai',
		sdf = "results/external-calls/evaluation/SDF"
	output:
		sample_vcf = "results/external-calls/evaluation/{truthset}/samples/{sample}-vcfeval.vcf.gz",
		fn = "results/external-calls/evaluation/{truthset}/samples/{sample}/vcfeval/fn.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	params:
		tmp = "results/external-calls/evaluation/{truthset}/samples/{sample}/vcfeval_temp",
		outname = "results/external-calls/evaluation/{truthset}/samples/{sample}/vcfeval"
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=59
	log:
		"results/external-calls/evaluation/{truthset}/samples/{sample}/vcfeval.log"
	shell:
		"""
		bcftools view --samples {wildcards.sample} {input.panel} | bcftools view --min-ac 1 | python3 workflow/scripts/fix-header.py {input.ref_index} | bgzip -c > {output.sample_vcf}
		tabix -p vcf {output.sample_vcf}
		rtg vcfeval -b {input.truthset} -c {output.sample_vcf} -t {input.sdf} -o {params.tmp} --squash-ploidy  &> {log}
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""



# same as previous rule, but for SVs using truvari instead of vcfeval
ne_false_negatives_truvari:
	input:
		truthset = "results/external-calls/evaluation/{truthset}/truthset-{truthset}_sv.vcf.gz",
		panel = PANEL_BI,
		reference = REFERENCE,
		ref_index = REFERENCE + '.fai',
	output:
		sample_vcf = "results/external-calls/evaluation/{truthset}/samples/{sample}-truvari.vcf.gz",
		fn = "results/external-calls/evaluation/{truthset}/samples/{sample}/truvari/fn.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	params:
		tmp = "results/external-calls/evaluation/{truthset}/samples/truvari_temp",
		outname = "results/external-calls/evaluation/{truthset}/samples/{sample}/truvari"
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=59
	log:
		"results/external-calls/evaluation/{truthset}/samples/{sample}/truvari.log"
	shell:
		"""
		bcftools view --samples {wildcards.sample} {input.panel} | bcftools view --min-ac 1 | python3 workflow/scripts/fix-header.py {input.ref_index} | bgzip -c > {output.sample_vcf}
		tabix -p vcf {output.sample_vcf}
		truvari bench -b {input.truthset} -c {output.sample_vcf} -f {input.reference} -o {params.tmp} --pick ac -r 2000 -C 5000 --passonly &> {log}
		bgzip {params.tmp}/fn.vcf
		tabix -p vcf {params.tmp}/fn.vcf.gz
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""


rule determine_unique:
	"""
	intersect the sets of FNs computed for each panel sample. The intersection then defines the set of unique/untypable variants
	"""
	input:
		expand("results/external-calls/evaluation/{{truthset}}/samples/{sample}/{{method}}/fn.vcf.gz", sample=PANEL_SAMPLES)
	output:
		
		unique_tsv = "results/external-calls/evaluation/{truthset}/{truthset}-unique_{method}.tsv",
		unique_vcf = "results/external-calls/evaluation/{{truthset}/{truthset}-unique_{method}.vcf"
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



