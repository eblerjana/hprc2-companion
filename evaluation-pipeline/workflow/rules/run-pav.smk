

rule prepare_assembly_table:
	input:
		cons1 = "{results}/evaluation/{callset}/{sample}_hap1/{sample}_hap1_consensus.fa.gz",
		cons2 = "{results}/evaluation/{callset}/{sample}_hap2/{sample}_hap2_consensus.fa.gz",
		assembly = lambda wildcards: ASSEMBLIES[wildcards.sample]['hap2'] if wildcards.haplotype == "hap1" else ASSEMBLIES[wildcards.sample]['hap1'],
		reference = lambda wildcards: ASSEMBLIES[wildcards.sample][wildcards.haplotype]

	output:
		tsv = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/assemblies.tsv",
		h1 = temp("{results}/evaluation/pav_{callset}_{sample}_{haplotype}/h1.fa.gz"),
		h2 = temp("{results}/evaluation/pav_{callset}_{sample}_{haplotype}/h2.fa.gz"),
		h3 = temp("{results}/evaluation/pav_{callset}_{sample}_{haplotype}/h3.fa.gz"),
		ref = temp("{results}/evaluation/pav_{callset}_{sample}_{haplotype}/ref.fa.gz"),
		fai = temp("{results}/evaluation/pav_{callset}_{sample}_{haplotype}/ref.fa.gz.fai")
	shell:
		"""
		cp {input.cons1} {output.h1}
		cp {input.cons2} {output.h2}
		cp {input.assembly} {output.h3}
		cp {input.reference} {output.ref}
		cp {input.reference}.fai {output.fai}
		echo \"NAME\tHAP_h1\tHAP_h2\tHAP_h3\" > {output.tsv}
		echo \"{wildcards.callset}_{wildcards.sample}_{wildcards.haplotype}\th1.fa.gz\th2.fa.gz\th3.fa.gz\" >> {output.tsv}
		"""


rule prepare_pav_contig:
	input:
		fasta = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/ref.fa.gz"
	output:
		json = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/config.json"
	run:
		with open(output.json, 'w') as outfile:
			outfile.write("{\n")
			outfile.write("\t\"reference\": \"ref.fa.gz\",\n")
			outfile.write("\t\"no_link_qry\": \"True\"\n")
			outfile.write("}\n")


rule run_pav:
	input:
		tsv = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/assemblies.tsv",
		pav_config = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/config.json",
		h1 = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/h1.fa.gz",
		h2 = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/h2.fa.gz",
		h3 = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/h3.fa.gz",
		ref = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/ref.fa.gz",
		fai = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/ref.fa.gz.fai"
	output:
		complete = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/run.complete",
		vcf = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/{callset}_{sample}_{haplotype}.vcf.gz"
	threads:
		24
	params:
		wdir = "{results}/evaluation/pav_{callset}_{sample}_{haplotype}/"
	log:
		"{results}/evaluation/pav_{callset}_{sample}_{haplotype}/pav.log"
	benchmark:
		"{results}/evaluation/pav_{callset}_{sample}_{haplotype}/pav.benchmark"
	resources:
		mem_mb = 150000,
		walltime = "24:00:00"	
	singularity:
		"workflow/container/pav_latest.sif"
	shell:
		"""
		snakemake --verbose --jobs {threads} -d {params.wdir} -s /opt/pav/Snakefile --rerun-incomplete --restart-times 0 --notemp &> {log} && touch {output.complete}
		"""
