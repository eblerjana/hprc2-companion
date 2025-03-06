

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
		temp("{results}/leave-one-out/preprocessed-vcfs/{sample}_{representation}_no-missing.vcf")
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


rule leave_one_out_prepare_panel_multi:
	"""
	Create input genotyping by removing evaluation sample from multi-allelic panel VCF.
	"""
	input:
		PANEL_MULTI
	output:
		temp("{results}/leave-one-out/input-panel/panel-{sample}_multi.vcf")
	conda:
		"../envs/genotyping.yml"
	priority: 1
	log:
		"{results}/leave-one-out/input-panel/panel-{sample}_multi.log"
	resources:
		mem_mb=20000
	shell:
		"""
		bcftools view --samples ^{wildcards.sample} {input} | bcftools view --min-ac 1 2> {log} 1> {output}
		"""


rule leave_one_out_prepare_panel_bi:
	"""
	Create input genotyping by removing evaluation sample from bi-allelic panel VCF.
	"""
	input:
		PANEL_BI
	output:
		temp("{results}/leave-one-out/input-panel/panel-{sample}_bi.vcf")
	conda:
		"../envs/genotyping.yml"
	priority: 1
	log:
		"{results}/leave-one-out/input-panel/panel-{sample}_bi.log"
	resources:
		mem_mb=20000
	shell:
		"""
		bcftools view --samples ^{wildcards.sample} {input} | bcftools view --min-ac 1 2> {log} 1> {output}
		"""

