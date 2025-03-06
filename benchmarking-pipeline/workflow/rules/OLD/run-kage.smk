
rule kage_index:
	input:
		reference = REFERENCE,
		VCF = PANEL_BI
	output:
		"{results}/kage/sampled-{size}/{sample}/{sample}-pangenie_multi_genotyping.vcf"	
	resources:
	benchmark:
	singularity:
	shell:
