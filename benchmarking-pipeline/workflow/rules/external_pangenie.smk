
rule external_pangenie_prepare_panel:
	input:
		vcf = PANEL_MULTI
	output:
		temp("{results}/external-calls/panel-multi.vcf")
	shell:
		"""
		gunzip -c {input} > {output}
		"""


########################################################
##################    run PanGenie    ##################
########################################################



rule external_pangenie_index:
	"""
	Create index for PanGenie.
	"""
	input:
		vcf = "{results}/external-calls/panel-multi.vcf",
		fasta = REFERENCE,
	output:
		temp(directory("{results}/external-calls/pangenie/index/"))
	log:
		"{results}/external-calls/pangenie/index.log"
	resources:
		mem_mb = 110000,
		walltime = "4:00:00"
	threads: 24
	params:
		out_prefix = "{results}/external-calls/pangenie/index/index"
	benchmark:
		"{results}/external-calls/pangenie/index.benchmark.txt"
	singularity:
		"workflow/container/pangenie.sif"
	shell:
		"""
		mkdir {output}
		PanGenie-index -v {input.vcf} -r {input.fasta} -o {params.out_prefix} -t {threads}  &> {log}
		"""


rule external_pangenie_genotype_subset:
	"""
	Run genotyping using the full panel (no sampling).
	"""
	input:
		reads = lambda wildcards: ILLUMINA[wildcards.sample],
		index = "{results}/external-calls/pangenie/index/"
	output:
		temp("{results}/external-calls/pangenie-subset-{size}/{sample}/{sample}_pangenie-subset-{size}_multi_genotyping.vcf")
	log:
		"{results}/external-calls/pangenie-subset-{size}/{sample}/{sample}_pangenie-subset-{size}_multi_genotyping.log"
	resources:
		mem_mb = 100000,
		walltime = "10:00:00"
	params:
		index = "{results}/external-calls/pangenie/index/index",
		out_prefix = "{results}/external-calls/pangenie-subset-{size}/{sample}/{sample}_pangenie-subset-{size}_multi"
	benchmark:
		"{results}/external-calls/pangenie-subset-{size}/{sample}/{sample}_pangenie-subset-{size}.benchmark.txt"
	threads:
		24
	singularity:
		"workflow/container/pangenie.sif"
	shell:
		"""
		PanGenie -f {params.index} -i <(gunzip -c {input.reads}) -o {params.out_prefix} -t {threads} -j {threads} -a {wildcards.size} &> {log}
		"""


rule external_pangenie_genotype_sampling:
	"""
	Run genotyping using sampling.
	"""
	input:
		reads = lambda wildcards: ILLUMINA[wildcards.sample],
		index = "{results}/external-calls/pangenie/index/"
	output:
		genotypes = temp("{results}/external-calls/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}_multi_genotyping.vcf"),
		panel = temp("{results}/external-calls/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}_multi_panel.vcf"),
		paths = temp(expand("{{results}}/external-calls/pangenie-sampled-{{size}}/{{sample}}/{{sample}}_pangenie-sampled-{{size}}_multi_paths_{chromosome}.tsv", chromosome = CHROMOSOMES))
	log:
		"{results}/external-calls/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}_multi_genotyping.log"
	resources:
		mem_mb = 60000,
		walltime = "3:00:00"
	wildcard_constraints:
		size = "[0-9,-.,x]+"
	params:
		index = "{results}/external-calls/pangenie/index/index",
		out_prefix = "{results}/external-calls/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}_multi",
		size = lambda wildcards: wildcards.size.split('-')[0],
		penalty = lambda wildcards: wildcards.size.split('-')[-1].split('x')[0],
		pop_size = lambda wildcards: wildcards.size.split('-')[-1].split('x')[1]
	benchmark:
		"{results}/external-calls/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}.benchmark.txt"
	threads:
		24
	singularity:
		"workflow/container/pangenie.sif"
	shell:
		"""
		PanGenie -f {params.index} -i <(gunzip -c {input.reads}) -o {params.out_prefix} -t {threads} -j {threads} -d -x {params.size} -y {params.penalty} -b {params.pop_size}  &> {log}
		"""




rule external_convert_genotypes_to_biallelic:
	"""
	Convert genotyped VCF to biallelic representation.
	"""
	input:
		vcf = "{results}/external-calls/pangenie-{mode}/{sample}/{sample}_pangenie-{mode}_multi_genotyping.vcf",
		biallelic = PANEL_BI
	output:
		bi = temp("{results}/external-calls/pangenie-{mode}/{sample}/{sample}_pangenie-{mode}_bi_genotyping.vcf.gz"),
		bi_tbi = temp("{results}/external-calls/pangenie-{mode}/{sample}/{sample}_pangenie-{mode}_bi_genotyping.vcf.gz.tbi"),
		multi = "{results}/external-calls/pangenie-{mode}/{sample}/{sample}_pangenie-{mode}_multi_genotyping.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb=30000
	priority: 1
	shell:
		"""
		bgzip -c {input.vcf} > {output.multi}
		zcat {output.multi} | python3 workflow/scripts/convert-to-biallelic.py {input.biallelic} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' | bgzip > {output.bi}
		tabix -p vcf {output.bi}
		"""


rule external_compress_pangenie_paths:
	input:
		"{results}/external-calls/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}_multi_paths_{chromosome}.tsv"
	output:
		"{results}/external-calls/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}_multi_paths_{chromosome}.tsv.gz"
	shell:
		"""
		bgzip -c {input} > {output}
		"""
