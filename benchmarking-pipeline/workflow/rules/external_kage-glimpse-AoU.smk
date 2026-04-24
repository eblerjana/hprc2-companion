

###############################################################################################################################################
######################################## postprocess externally produced KAGE-Glimopse results  ###############################################
###############################################################################################################################################


rule external_kage_AoU_bubble:
	"""
	Prepare Kage VCF in bubble representation
	"""
	input:
		lambda wildcards: KAGE_AOU if wildcards.genotyper == "kage-AoU-bubble" else KAGE_AOU_POPPED
	output:
		"{results}/external-calls/{genotyper}/{sample}/{sample}_{genotyper}_bi_genotyping.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		genotyper = "kage-AoU-bubble|kage-AoU-popped"
	resources:
		mem_mb = 50000
	shell:
		"""
		bcftools view -s {wildcards.sample} {input} -Oz -o {output}
		tabix -p vcf {output}
		"""


################################### annotate the variants and decompose bubbles based on annotations  #########################################


rule external_kage_AoU_prepare_panel:
	"""
	Prepare biallelic representation of multiallelic
	panel.
	"""
	input:
		PANEL_MULTI
	output:
		"{results}/external-calls/kage-AoU/panel.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb = 50000,
		walltime = "05:00:00"
	shell:
		"""
		bcftools norm -m- {input} | bgzip > {output}
		tabix -p vcf {output}
		"""


rule external_kage_AoU_postprocess:
	"""
	Kage does not keep track of the variant IDs in the INFO field.
	This rule adds them back in (assuming the variant representation is identical)
	"""
	input:
		kage = KAGE_AOU,
		vcf = "{results}/external-calls/kage-AoU/panel.vcf.gz"
	output:
		vcf = temp("{results}/external-calls/kage-AoU/{sample}/{sample}_kage-AoU_bi_genotyping_tmp.vcf.gz"),
		tbi = "{results}/external-calls/kage-AoU/{sample}/{sample}_kage-AoU_bi_genotyping_tmp.vcf.gz.tbi"
	resources:
		mem_mb = 100000,
		walltime = "15:00:00"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools view -s {wildcards.sample} {input.kage} | python3 workflow/scripts/add-ids.py {input.vcf} {wildcards.sample} | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


rule external_kage_AoU_convert_to_bi:
	"""
	Convert VCF to biallelic representation based on
	annotations added by previous rule.
	"""
	input:
		kage = "{results}/external-calls/kage-AoU/{sample}/{sample}_kage-AoU_bi_genotyping_tmp.vcf.gz",
		panel = PANEL_BI
	output:
		 "{results}/external-calls/kage-AoU/{sample}/{sample}_kage-AoU_bi_genotyping.vcf.gz"
	resources:
		mem_mb = 50000
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		zcat {input.kage} | python3 workflow/scripts/convert-to-biallelic.py {input.panel} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' | bgzip > {output}
		tabix -p vcf {output}
		"""
