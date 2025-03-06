

################################################################
######   prepare input panel and ground truth genotypes  #######
################################################################


rule leave_one_out_remove_missing:
	"""
	Remove positions that are ".|." in the left out sample. These cannot be used for evaluation, as the true genotype is unknown
	for now, also remove CHM13, because PanGenie cannot handle haploid samples
	remove sample from the panel.
	"""
	input:
		lambda wildcards: PANEL_MULTI if wildcards.representation == 'multi' else PANEL_BI
	output:
		temp("{results}/leave-one-out/pangenie/preprocessed-vcfs/{sample}_{representation}_no-missing.vcf")
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb=20000
	priority: 1
	wildcard_constraints:
		representation = "bi|multi"
	shell:
		"""
		zcat {input} | python3 workflow/scripts/remove-missing.py {wildcards.sample} > {output}
		"""


rule leave_one_out_prepare_panel:
	"""
	Create input genotyping by removing evaluation sample from panel VCF.
	"""
	input:
		PANEL_MULTI
	output:
		temp("{results}/leave-one-out/pangenie/input-panel/panel-{sample}.vcf")
	conda:
		"../envs/genotyping.yml"
	priority: 1
	log:
		"{results}/leave-one-out/pangenie/input-panel/panel-{sample}.log"
	resources:
		mem_mb=20000
	shell:
		"""
		bcftools view --samples ^{wildcards.sample} {input} | bcftools view --min-ac 1 2> {log} 1> {output}
		"""



########################################################
##################    run PanGenie    ##################
########################################################



rule leave_one_out_pangenie_index:
	"""
	Create index for PanGenie.
	"""
	input:
		vcf = "{results}/leave-one-out/pangenie/input-panel/panel-{sample}.vcf",
		fasta = REFERENCE,
	output:
		temp(directory("{results}/leave-one-out/pangenie/index-{sample}/"))
	log:
		"{results}/leave-one-out/pangenie/index-{sample}.log"
	resources:
		mem_mb = 100000,
		walltime = "5:00:00"
	threads: 24
	params:
		out_prefix = "{results}/leave-one-out/pangenie/index-{sample}/index"
	benchmark:
		"{results}/leave-one-out/pangenie/index-{sample}/index-{sample}_benchmark.txt"
	singularity:
		"workflow/container/eblerjana_eblerjana_pangenie-sampler.sif"
	shell:
		"""
		PanGenie-index -v {input.vcf} -r {input.fasta} -o {params.out_prefix} -t {threads}  &> {log}
		"""


rule leave_one_out_pangenie_genotype_subset:
	"""
	Run genotyping using the full panel (no sampling).
	"""
	input:
		reads = lambda wildcards: READS[wildcards.sample],
		index = "{results}/leave-one-out/pangenie/index-{sample}/"
	output:
		temp("{results}/leave-one-out/pangenie/subset-{size}/{sample}/{sample}-pangenie_multi_genotyping.vcf")
	log:
		"{results}/leave-one-out/pangenie/subset-{size}/{sample}/{sample}-pangenie_multi_genotyping.log"
	resources:
		mem_mb = 100000,
		walltime = "10:00:00"
	params:
		index = "{results}/leave-one-out/pangenie/index-{sample}/index",
		out_prefix = "{results}/leave-one-out/pangenie/subset-{size}/{sample}/{sample}-pangenie_multi"
	benchmark:
		"{results}/leave-one-out/pangenie/subset-{size}/{sample}/{sample}-pangenie_multi_benchmark.txt"
	threads:
		24
	singularity:
		"workflow/container/eblerjana_eblerjana_pangenie-sampler.sif"
	shell:
		"""
		PanGenie -f {params.index} -i <(gunzip -c {input.reads}) -o {params.out_prefix} -t {threads} -j {threads} -a {wildcards.size} &> {log}
		"""


rule leave_one_out_pangenie_genotype_sampling:
	"""
	Run genotyping using sampling.
	"""
	input:
		reads = lambda wildcards: READS[wildcards.sample],
		index = "{results}/leave-one-out/pangenie/index-{sample}/"
	output:
		genotypes = temp("{results}/leave-one-out/pangenie/sampled-{size}/{sample}/{sample}-pangenie_multi_genotyping.vcf"),
		panel = "{results}/leave-one-out/pangenie/sampled-{size}/{sample}/{sample}-pangenie_multi_panel.vcf",
		paths = expand("{{results}}/leave-one-out/pangenie/sampled-{{size}}/{{sample}}/{{sample}}-pangenie_multi_paths_{chromosome}.tsv", chromosome = CHROMOSOMES)
	log:
		"{results}/leave-one-out/pangenie/sampled-{size}/{sample}/{sample}-pangenie_multi_genotyping.log"
	resources:
		mem_mb = 100000,
		walltime = "5:00:00"
	wildcard_constraints:
		size = "[0-9,-.,x]+"
	params:
		index = "{results}/leave-one-out/pangenie/index-{sample}/index",
		out_prefix = "{results}/leave-one-out/pangenie/sampled-{size}/{sample}/{sample}-pangenie_multi",
		size = lambda wildcards: wildcards.size.split('-')[0],
		penalty = lambda wildcards: wildcards.size.split('-')[-1].split('x')[0],
		pop_size = lambda wildcards: wildcards.size.split('-')[-1].split('x')[1]
	benchmark:
		"{results}/leave-one-out/pangenie/sampled-{size}/{sample}/{sample}-pangenie_multi_benchmark.txt"
	threads:
		24
	singularity:
		"workflow/container/eblerjana_eblerjana_pangenie-sampler.sif"
	shell:
		"""
		PanGenie -f {params.index} -i <(gunzip -c {input.reads}) -o {params.out_prefix} -t {threads} -j {threads} -d -x {params.size} -y {params.penalty} -b {params.pop_size}  &> {log}
		"""




rule leave_one_out_convert_genotypes_to_biallelic:
	"""
	Convert genotyped VCF to biallelic representation.
	"""
	input:
		vcf = "{results}/leave-one-out/pangenie/{mode}/{sample}/{sample}-pangenie_multi_genotyping.vcf",
		biallelic = PANEL_BI
	output:
		bi = temp("{results}/leave-one-out/pangenie/{mode}/{sample}/{sample}-pangenie_bi_genotyping.vcf"),
		multi = "{results}/leave-one-out/pangenie/{mode}/{sample}/{sample}-pangenie_multi_genotyping.vcf.gz"
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


