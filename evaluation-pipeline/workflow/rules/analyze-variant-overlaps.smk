BUBBLE_SIZE_THRESHOLD = 50

rule overlaps_copy_reference_data:
	"""
	Create local copies of data for PAV
	"""
	input:
		assemb1 = lambda wildcards: ASSEMBLIES[wildcards.sample]['hap1'],
		assemb2 = lambda wildcards: ASSEMBLIES[wildcards.sample]['hap2'],
		reference = lambda wildcards: OVERLAPS[wildcards.set]["reference"]
	output:
		assemb1 = temp("{results}/variant-overlaps/{set}/pav_{set}_{sample}/assembly_hap1.fa.gz"),
		assemb2 = temp("{results}/variant-overlaps/{set}/pav_{set}_{sample}/assembly_hap2.fa.gz"),
		reference = temp("{results}/variant-overlaps/{set}/pav_{set}_{sample}/ref.fa"),
		fai = temp("{results}/variant-overlaps/{set}/pav_{set}_{sample}/ref.fa.fai")
	shell:
		"""
		cp {input.assemb1} {output.assemb1}
		cp {input.assemb2} {output.assemb2}
		cp {input.reference} {output.reference}
		cp {input.reference}.fai {output.fai}
		"""

rule overlaps_copy_haplotype_data:
	"""
	Create local copy of consensus haplotypes for PAV
	"""
	input:
		"{results}/evaluation/{callset}/{sample}_{haplotype}/{sample}_{haplotype}_consensus.fa.gz"
	output:
		temp("{results}/variant-overlaps/{set}/pav_{set}_{sample}/{callset}_{haplotype}.fa.gz")
	shell:
		"""
		cp {input} {output}
		"""

rule overlaps_assembly_table:
	"""
	Create input file for PAV
	"""
	output:
		tsv = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/assemblies.tsv"
	params:
		callsets = lambda wildcards: OVERLAPS[wildcards.set]["callsets"]
	run:
		with open(output.tsv, 'w') as outfile:
			header = ["NAME"]
			for i in range(0, len(params.callsets) * 2 + 2):
				header.append("HAP_h" + str(i+1))
			outfile.write("\t".join(header) + "\n")
			filenames = [wildcards.set, "assembly_hap1.fa.gz", "assembly_hap2.fa.gz"]
			for callset in params.callsets:
				filenames.append(callset + "_hap1.fa.gz")
				filenames.append(callset + "_hap2.fa.gz")
			outfile.write("\t".join(filenames) + "\n")
		
	
rule overlaps_prepare_pav_contig:
	input:
		"{results}/variant-overlaps/{set}/pav_{set}_{sample}/ref.fa"
	output:
		json = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/config.json"
	run:
		with open(output.json, 'w') as outfile:
			outfile.write("{\n")
			outfile.write("\t\"reference\": \"ref.fa\",\n")
			outfile.write("\t\"no_link_qry\": \"True\"\n")
			outfile.write("}\n")


rule overlaps_run_pav:
	input:
		tsv = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/assemblies.tsv",
		pav_config = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/config.json",
		assemb1 = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/assembly_hap1.fa.gz",
		assemb2 = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/assembly_hap2.fa.gz",
		reference = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/ref.fa",
		fai = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/ref.fa.fai",
		consensus = lambda wildcards: expand("{{results}}/variant-overlaps/{{set}}/pav_{{set}}_{{sample}}/{callset}_{haplotype}.fa.gz", callset = OVERLAPS[wildcards.set]["callsets"], haplotype = ["hap1", "hap2"] )
	output:
		complete = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/run.complete",
		vcf = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/{set}.vcf.gz"
	threads:
		24
	params:
		wdir = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/"
	log:
		"{results}/variant-overlaps/{set}/pav_{set}_{sample}/pav.log"
	benchmark:
		"{results}/variant-overlaps/{set}/pav_{set}_{sample}/pav.benchmark"
	resources:
		mem_mb = 150000,
		walltime = "24:00:00"	
	singularity:
		"workflow/container/pav_latest.sif"
	shell:
		"""
		snakemake --verbose --jobs {threads} -d {params.wdir} -s /opt/pav/Snakefile --rerun-incomplete --restart-times 0 --notemp &> {log} && touch {output.complete}
		"""


rule overlaps_filter_vcf:
	"""
	Only keep variants with filter PASS.
	"""
	input:
		vcf = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/{set}.vcf.gz",
		reference = lambda wildcards: OVERLAPS[wildcards.set]["reference"]
	output:
		 "{results}/variant-overlaps/{set}/pav_{set}_{sample}/{set}_pass.vcf.gz"
	conda:
		"../envs/overlaps.yml"
	resources:
		mem_mb = 50000
	shell:
		"""
		bcftools view -f PASS {input.vcf} | bcftools norm -m -any -f {input.reference} | bcftools sort -Oz -o {output}
		tabix -p vcf {output} 
		"""


rule overlaps_filter_callable_vcf:
	"""
	Only keep variants inside of callable regions provided.
	"""
	input:
		vcf = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/{set}_pass.vcf.gz",
		bed = lambda wildcards: OVERLAPS[wildcards.set]["callable_regions"],
		reference = lambda wildcards: OVERLAPS[wildcards.set]["reference"]
	output:
		 "{results}/variant-overlaps/{set}/pav_{set}_{sample}/{set}_callable.vcf.gz"
	conda:
		"../envs/overlaps.yml"
	shell:
		"""
		bedtools subtract -a {input.vcf} -b {input.bed} -header | bgzip > {output}
		tabix -p vcf {output} 
		"""


rule overlaps_get_bubble_coordinates:
	"""
	Write a BED file with coordinates of all bubbles of
	at least BUBBLE_SIZE_THRESHOLD.
	"""
	input:
		lambda wildcards: OVERLAPS[wildcards.set]["bubbles"]
	output:
		"{results}/variant-overlaps/{set}/bubbles.bed"
	shell:
		"""
		zcat {input} | python3 workflow/scripts/vcf-to-bed.py {BUBBLE_SIZE_THRESHOLD} > {output}
		"""


rule overlaps_plot_overlaps:
	"""
	Analyze and plot overlaps of variants across compared
	callsets.
	"""
	input:
		vcf = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/{set}_{filter}.vcf.gz",
		bed =  "{results}/variant-overlaps/{set}/bubbles.bed",
		pav_config = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/assemblies.tsv",
		variants = "{results}/variant-overlaps/{set}/evaluation/truvari/truvari_{set}_{sample}_{filter}/tp-comp.vcf.gz",
		closest = "{results}/variant-overlaps/{set}/evaluation/truvari/{set}_{sample}_{filter}_closest.tsv"
	output:
		tsv = "{results}/variant-overlaps/{set}/evaluation/{set}_{sample}_{filter}.tsv",
		pdf = "{results}/variant-overlaps/{set}/evaluation/{set}_{sample}_{filter}.pdf"
	conda:
		"../envs/overlaps.yml"
	resources:
		mem_mb = 50000
	log:
		"{results}/variant-overlaps/{set}/evaluation/{set}_{sample}_{filter}.log"
	params:
		outname = "{results}/variant-overlaps/{set}/evaluation/{set}_{sample}_{filter}"
	shell:
		"""
		bedtools annotate -i {input.vcf} -files {input.bed} | python3 workflow/scripts/count-variants.py --config {input.pav_config} --outname {params.outname} -t {BUBBLE_SIZE_THRESHOLD} -v {input.variants} -n {input.closest} &> {log}
		"""


######### determine variant overlaps ########

rule overlaps_truvari_prepare_bubbles:
	"""
	Normalize graph variants, set their genotypes to 1/1 (to make truvari
	consider all of them), extract SVs and fix format to work with
	truvari.
	"""
	input:
		vcf = lambda wildcards: OVERLAPS[wildcards.set]["variants"],
		reference = lambda wildcards: OVERLAPS[wildcards.set]["reference"],
		ref_index = lambda wildcards: OVERLAPS[wildcards.set]["reference"] + ".fai"
	output:
		"{results}/variant-overlaps/{set}/evaluation/truvari/{set}_graph.vcf.gz"
	conda:
		"../envs/truvari.yml"
	resources:
		mem_mb = 50000,
		walltime = "05:00:00"
	shell:
		"""
		bcftools norm -m -any -f {input.reference} {input.vcf} | python3 workflow/scripts/set-hom-gt.py | python3 workflow/scripts/fix-header.py {input.ref_index} | python3 workflow/scripts/extract-varianttype.py sv | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' | bgzip > {output}
		tabix -p vcf {output}
		"""

	
rule overlaps_truvari_prepare_calls:
	"""
	Normalize calls, set their genotypes to 1/1, remove INV calls, extract SVs
	and fix format to work with truvari.
	"""
	input:
		vcf = "{results}/variant-overlaps/{set}/pav_{set}_{sample}/{set}_{filter}.vcf.gz",
		reference = lambda wildcards: OVERLAPS[wildcards.set]["reference"],
		ref_index = lambda wildcards: OVERLAPS[wildcards.set]["reference"] + ".fai"
	output:
		"{results}/variant-overlaps/{set}/evaluation/truvari/{set}_{sample}_{filter}_calls.vcf.gz"
	conda:
		"../envs/truvari.yml"
	resources:
		mem_mb = 10000
	shell:
		"""
		zcat {input.vcf} | python3 workflow/scripts/set-hom-gt.py | python3 workflow/scripts/fix-header.py {input.ref_index} | grep -v "\<INV\>" | python3 workflow/scripts/extract-varianttype.py sv | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' | bgzip > {output}
		tabix -p vcf {output}
		"""


rule overlaps_truvari:
	input:
		callset = "{results}/variant-overlaps/{set}/evaluation/truvari/{set}_{sample}_{filter}_calls.vcf.gz",
		baseline = "{results}/variant-overlaps/{set}/evaluation/truvari/{set}_graph.vcf.gz",
		reference = lambda wildcards: OVERLAPS[wildcards.set]["reference"]
	output:
		tp = "{results}/variant-overlaps/{set}/evaluation/truvari/truvari_{set}_{sample}_{filter}/tp-comp.vcf.gz"	
	conda:
		"../envs/truvari.yml"
	log:
		"{results}/variant-overlaps/{set}/evaluation/truvari/truvari_{set}_{sample}_{filter}.log"
	params:
		tmp = "{results}/variant-overlaps/{set}/evaluation/truvari/truvari_{set}_{sample}_{filter}_tmp",
		outname = "{results}/variant-overlaps/{set}/evaluation/truvari/truvari_{set}_{sample}_{filter}"
	resources:
		mem_mb=20000,
		walltime = "00:40:00"
	shell:
		"""
		truvari bench -b {input.baseline} -c {input.callset} -f {input.reference} -o {params.tmp} --pick multi -r 1000 -C 1000 -s 50 -S 15 --sizemax 100000 -p 0.0 -P 0.3 -O 0.0 --passonly --no-ref a &> {log}
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""

rule overlaps_closest:
	"""
	For all SV calls, compute distance to the closest graph SV.
	"""
	input:
		calls = "{results}/variant-overlaps/{set}/evaluation/truvari/{set}_{sample}_{filter}_calls.vcf.gz",
		graph = "{results}/variant-overlaps/{set}/evaluation/truvari/{set}_graph.vcf.gz"
	output:
		"{results}/variant-overlaps/{set}/evaluation/truvari/{set}_{sample}_{filter}_closest.tsv"
	conda:
		 "../envs/overlaps.yml"
	resources:
		mem_mb=20000
	shell:
		"""
		bedtools closest -a {input.calls} -b {input.graph} -d -t first > {output}
		"""
