
########################################################
##################    run PanGenie    ##################
########################################################



rule leave_one_out_pangenie_index:
	"""
	Create index for PanGenie.
	"""
	input:
		vcf = "{results}/leave-one-out/input-panel/panel-{sample}_multi.vcf",
		fasta = REFERENCE,
	output:
		temp(directory("{results}/leave-one-out/pangenie/index-{sample}/"))
	log:
		"{results}/leave-one-out/pangenie/index-{sample}.log"
	resources:
		mem_mb = 110000,
		walltime = "4:00:00"
	threads: 24
	params:
		out_prefix = "{results}/leave-one-out/pangenie/index-{sample}/index"
	benchmark:
		"{results}/leave-one-out/pangenie/index-{sample}.benchmark.txt"
	singularity:
		"workflow/container/pangenie.sif"
	shell:
		"""
		mkdir {output}
		PanGenie-index -v {input.vcf} -r {input.fasta} -o {params.out_prefix} -t {threads}  &> {log}
		"""


rule leave_one_out_pangenie_genotype_subset:
	"""
	Run genotyping using the full panel (no sampling).
	"""
	input:
		reads = lambda wildcards: ILLUMINA[wildcards.sample],
		index = "{results}/leave-one-out/pangenie/index-{sample}/"
	output:
		temp("{results}/leave-one-out/pangenie-subset-{size}/{sample}/{sample}_pangenie-subset-{size}_multi_genotyping.vcf")
	log:
		"{results}/leave-one-out/pangenie-subset-{size}/{sample}/{sample}_pangenie-subset-{size}_multi_genotyping.log"
	resources:
		mem_mb = 100000,
		walltime = "10:00:00"
	params:
		index = "{results}/leave-one-out/pangenie/index-{sample}/index",
		out_prefix = "{results}/leave-one-out/pangenie-subset-{size}/{sample}/{sample}_pangenie-subset-{size}_multi"
	benchmark:
		"{results}/leave-one-out/pangenie-subset-{size}/{sample}/{sample}_pangenie-subset-{size}.benchmark.txt"
	threads:
		24
	singularity:
		"workflow/container/pangenie.sif"
	shell:
		"""
		PanGenie -f {params.index} -i <(gunzip -c {input.reads}) -o {params.out_prefix} -t {threads} -j {threads} -a {wildcards.size} &> {log}
		"""


rule leave_one_out_pangenie_genotype_sampling:
	"""
	Run genotyping using sampling.
	"""
	input:
		reads = lambda wildcards: ILLUMINA[wildcards.sample],
		index = "{results}/leave-one-out/pangenie/index-{sample}/"
	output:
		genotypes = temp("{results}/leave-one-out/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}_multi_genotyping.vcf"),
		panel = temp("{results}/leave-one-out/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}_multi_panel.vcf"),
		paths = temp(expand("{{results}}/leave-one-out/pangenie-sampled-{{size}}/{{sample}}/{{sample}}_pangenie-sampled-{{size}}_multi_paths_{chromosome}.tsv", chromosome = CHROMOSOMES))
	log:
		"{results}/leave-one-out/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}_multi_genotyping.log"
	resources:
		mem_mb = 60000,
		walltime = "3:00:00"
	wildcard_constraints:
		size = "[0-9,-.,x]+"
	params:
		index = "{results}/leave-one-out/pangenie/index-{sample}/index",
		out_prefix = "{results}/leave-one-out/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}_multi",
		size = lambda wildcards: wildcards.size.split('-')[0],
		penalty = lambda wildcards: wildcards.size.split('-')[-1].split('x')[0],
		pop_size = lambda wildcards: wildcards.size.split('-')[-1].split('x')[1]
	benchmark:
		"{results}/leave-one-out/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}.benchmark.txt"
	threads:
		24
	singularity:
		"workflow/container/pangenie.sif"
	shell:
		"""
		PanGenie -f {params.index} -i <(gunzip -c {input.reads}) -o {params.out_prefix} -t {threads} -j {threads} -d -x {params.size} -y {params.penalty} -b {params.pop_size}  &> {log}
		"""




rule leave_one_out_convert_genotypes_to_biallelic:
	"""
	Convert genotyped VCF to biallelic representation.
	"""
	input:
		vcf = "{results}/leave-one-out/pangenie-{mode}/{sample}/{sample}_pangenie-{mode}_multi_genotyping.vcf",
		biallelic = PANEL_BI
	output:
		bi = temp("{results}/leave-one-out/pangenie-{mode}/{sample}/{sample}_pangenie-{mode}_bi_genotyping.vcf"),
		multi = "{results}/leave-one-out/pangenie-{mode}/{sample}/{sample}_pangenie-{mode}_multi_genotyping.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb=30000
	priority: 1
	shell:
		"""
		bgzip -c {input.vcf} > {output.multi}
		zcat {output.multi} | python3 workflow/scripts/convert-to-biallelic.py {input.biallelic} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output.bi}
		"""


rule leave_one_out_compress_pangenie_paths:
	input:
		"{results}/leave-one-out/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}_multi_paths_{chromosome}.tsv"
	output:
		"{results}/leave-one-out/pangenie-sampled-{size}/{sample}/{sample}_pangenie-sampled-{size}_multi_paths_{chromosome}.tsv.gz"
	shell:
		"""
		bgzip -c {input} > {output}
		"""
