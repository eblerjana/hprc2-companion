

###############################################################################################################################################
################################################################## run KAGE ###################################################################
###############################################################################################################################################



rule external_kage_prepare_panel:
	"""
	Prepare VCF by adding AC,AF,AN annotations.
	"""
	input:
		PANEL_BI
	output:
		vcf = "{results}/external-calls/panel-bi.vcf",
		gz = "{results}/external-calls/panel-bi.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb = 20000
	shell:
		"""
		bcftools +fill-tags {input} -Oz -o {output.gz} -- -t AC,AN,AF
		tabix -p vcf {output.gz}
		gunzip -c {output.gz} > {output.vcf}
		"""


rule external_kage_prepare_fasta:
	"""
	KAGE cannot handle extra strings after the chromosome name,
	even when separated by a whitespace. This rule removes such
	strings.
	"""
	input:
		REFERENCE
	output:
		"{results}/external-calls/kage/reference.fa"
	benchmark:
		"{results}/external-calls/kage/reference.benchmark"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		cut -d" " -f1 {input} > {output}
		samtools faidx {output}
		"""


rule external_kage_extract_chromosome_from_vcf:
	"""
        Extract a chromosome from the VCF.
	"""
	input:
		"{results}/external-calls/panel-bi.vcf.gz"
	output:
		vcf = temp("{results}/external-calls/panel-bi_{chrom}.vcf"),
		gz = "{results}/external-calls/panel-bi_{chrom}.vcf.gz"
	benchmark:
		"{results}/external-calls/panel-bi_{chrom}.benchmark"
	resources:
		mem_mb = 50000
	wildcard_constraints:
		chrom = "|".join(KAGE_CHROMOSOMES)
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools view -r {wildcards.chrom} {input} > {output.vcf}
		bgzip -c {output.vcf} > {output.gz}
		tabix -p vcf {output.gz}
		"""


rule external_kage_index:
	"""
	Run Kage indexing.
	"""
	input:
		reference = "{results}/external-calls/kage/reference.fa",
		vcf = "{results}/external-calls/panel-bi_{chrom}.vcf"
	output:
		"{results}/external-calls/kage/index/index_{chrom}.npz"
	log:
		"{results}/external-calls/kage/index/index_{chrom}.log"
	resources:
		mem_mb = 500000,
		walltime = "06:00:00"
	benchmark:
		"{results}/external-calls/kage/index/index_{chrom}.benchmark.txt"
	singularity:
		"workflow/container/kage.sif"
	wildcard_constraints:
		chrom = "|".join(KAGE_CHROMOSOMES)
	threads:
		24
	shell:
		"""
		kage index -r {input.reference} -v {input.vcf} -o {output} -k 31 -t {threads}  &> {log}
		"""


rule external_kage_genotype:
	"""
	Run Kage genotyping.
	"""
	input:
		panel = "{results}/external-calls/panel-bi_{chrom}.vcf.gz",
		reads = lambda wildcards: ILLUMINA[wildcards.sample],
		index = "{results}/external-calls/kage/index/index_{chrom}.npz"
	output:
		"{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping_tmp_{chrom}.vcf"
	log:
		"{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping_{chrom}.log"
	resources:
		mem_mb = 200000,
		walltime = "05:00:00"
	benchmark:
		 "{results}/external-calls/kage/{sample}/{sample}_kage_{chrom}.benchmark.txt"
	singularity:
		"workflow/container/kage.sif"
	wildcard_constraints:
		chrom = "|".join(KAGE_CHROMOSOMES)
	threads: 24
	shell:
		"""
		kage genotype -i {input.index} -r {input.reads} -t {threads} --average-coverage {AVG_ILLUMINA_COV} -k 31  -o {output} &> {log}
		"""


rule external_kage_postprocess:
	"""
	Kage does not keep track of the variant IDs in the INFO field.
	This rule adds them back in (assuming the variant representation is identical)
	"""
	input:
		kage = "{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping_tmp_{chrom}.vcf",
		vcf = "{results}/external-calls/panel-bi.vcf.gz"
	output:
		vcf = "{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping_{chrom}.vcf.gz",
		tbi = "{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping_{chrom}.vcf.gz.tbi"
	resources:
		mem_mb = 50000,
		walltime = "05:00:00"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		cat {input.kage} | python3 workflow/scripts/add-ids.py {input.vcf} {wildcards.sample} | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""





###############################################################################################################################################
###################################################### run Glimpse on the KAGE genotypes ######################################################
###############################################################################################################################################

rule external_kage_glimpse_prepare_genotypes:
	"""
	Prepare KAGE genotypes. This rule removes records for which
	at least one GL entry is -inf. These cause errors in Glimpse
	when running with the --input-GL field.
	"""
	input:
		genotyped_vcf = "{results}/external-calls/kage/{sample}/{sample}_kage_bi_genotyping_{chrom}.vcf.gz"
	output:
		"{results}/external-calls/kage/{sample}/{sample}_kage_noINF_{chrom}.vcf.gz"
	conda:
		"../envs/glimpse.yml"
	shell:
		"""
		zcat {input.genotyped_vcf} | sed -e 's/nan/-1000000.0/g' | sed -e 's/-inf/-1000000.0/g' | sed -e 's/inf/-1000000.0/g' | bgzip > {output}
		tabix -p vcf {output} 
		"""



rule external_kage_glimpse_prepare_panel:
	"""
	Prepare panel for Glimpse. This rule removes
	variants with missing genotype information.
	"""
	input:
		PANEL_BI
	output:
		"{results}/external-calls/kage-glimpse/panel/glimpse-bi_{chrom}.vcf.gz"
	conda:
		"../envs/glimpse.yml"
	resources:
		mem_mb = 50000
	shell:
		"""
		bcftools view -r {wildcards.chrom} {input} | bcftools +fill-tags -Oz -o {output} -- -t AN,AC
		tabix -p vcf {output}
		"""


rule external_kage_glimpse_chunk:
	"""
	Split the panel into smaller chunks.
	"""
	input:
		vcf = "{results}/external-calls/kage-glimpse/panel/glimpse-bi_{chrom}.vcf.gz"
	output:
		"{results}/external-calls/kage-glimpse/chunks/chunks_{chrom}.txt"
	log:
		"{results}/external-calls/kage-glimpse/chunks/chunks_{chrom}.log"
	benchmark:
		"{results}/external-calls/kage-glimpse/chunks/chunks_{chrom}.benchmark"
	conda:
		"../envs/glimpse.yml"
	threads: 24
	shell:
		"""
		GLIMPSE_chunk --input {input.vcf} --region {wildcards.chrom} --window-size 1000000 --window-count 1000 --buffer-size 250000 --buffer-count 250 --thread {threads} --output {output} &> {log}
		"""


checkpoint external_kage_glimpse_prepare_regions:
	"""
	Prepare chunk regions.
	"""
	input:
		"{results}/external-calls/kage-glimpse/chunks/chunks_{chrom}.txt"
	output:
		directory("{results}/external-calls/kage-glimpse/{sample}/split_{chrom}/")
	benchmark:
		"{results}/external-calls/kage-glimpse/{sample}/split_{chrom}.bechmark"
	conda:
		"../envs/glimpse.yml"
	wildcard_constraints:
		chrom = "|".join(KAGE_CHROMOSOMES)
	threads: 1
	resources:
		mem_mb = "10000"
	shell:
		"""
		mkdir -p {output}
		while IFS=\"\" read -r LINE || [ -n \"$LINE\" ];
		do
			printf -v ID \"%02d\" $(echo $LINE | cut -d\" \" -f1)
			IRG=$(echo $LINE | cut -d\" \" -f3)
			ORG=$(echo $LINE | cut -d\" \" -f4)
			echo ${{IRG}} ${{ORG}} > {output}/split_{wildcards.chrom}_${{IRG}}_${{ORG}}.txt
		done < {input}
		"""


rule external_kage_glimpse_phase:
	"""
	Phase variants using Glimpse.
	"""
	input:
#		map = lambda wildcards: MAPS[wildcards.chrom],
		reference_panel = "{results}/external-calls/kage-glimpse/panel/glimpse-bi_{chrom}.vcf.gz",
		genotyped_vcf = "{results}/external-calls/kage/{sample}/{sample}_kage_noINF_{chrom}.vcf.gz"
	output:
		temp("{results}/external-calls/kage-glimpse/{sample}/{sample}_kage-glimpse_bi_genotyping_{chrom}_{irg}_{org}.bcf")
	benchmark:
		"{results}/external-calls/kage-glimpse/{sample}/{sample}_kage-glimpse_bi_genotyping_{chrom}_{irg}_{org}.benchmark"
	log:
		"{results}/external-calls/kage-glimpse/{sample}/{sample}_kage-glimpse_bi_genotyping_{chrom}_{irg}_{org}.log"
	conda:
		"../envs/glimpse.yml"
	wildcard_constraints:
		chrom = "|".join(KAGE_CHROMOSOMES)
	params:
		map = lambda wildcards: "--map " +  MAPS[wildcards.chrom] if wildcards.chrom in MAPS else " "
	resources:
		mem_mb = 10000,
		walltime = "00:10:00"
	threads: 1
	shell:
		"""
		GLIMPSE_phase --input {input.genotyped_vcf} --input-GL --thread {threads} {params.map} --input-region {wildcards.irg}  --output-region {wildcards.org} --reference {input.reference_panel} --output {output} &> {log}
		"""


def aggregate_concat_phased_variants(wildcards):
        checkpoint_output = checkpoints.external_kage_glimpse_prepare_regions.get(**wildcards).output[0]
        return expand("{results}/external-calls/kage-glimpse/{sample}/{sample}_kage-glimpse_bi_genotyping_{chrom}_{regions}.bcf",
						results = wildcards.results,
						sample = wildcards.sample,
						chrom = wildcards.chrom,
						regions = glob_wildcards(os.path.join(checkpoint_output, "split_" + wildcards.chrom + "_{regions}.txt")).regions)


rule external_kage_glimpse_ligate:
	"""
	Combine phased chunks into full chromosome VCF.
	"""
	input:
		aggregate_concat_phased_variants
	output:
		lst = "{results}/external-calls/kage-glimpse/{sample}/glimpse_{chrom}.lst",
		vcf = temp("{results}/external-calls/kage-glimpse/{sample}/{sample}_kage-glimpse_bi_genotyping_tmp_{chrom}.vcf.gz")
	benchmark:
		"{results}/external-calls/kage-glimpse/{sample}/{sample}_kage-glimpse_bi_genotyping_{chrom}.benchmark"
	log:
		"{results}/external-calls/kage-glimpse/{sample}/{sample}_kage-glimpse_bi_genotyping_{chrom}.log"
	conda:
		"../envs/glimpse.yml"
	threads: 24
	resources:
		mem_mb = 10000
	shell:
		"""
		ls -1v {input} > {output.lst}
		GLIMPSE_ligate --input {output.lst} --thread {threads} --output {output.vcf} &> {log}
		tabix -p vcf {output.vcf}
		"""


rule external_kage_glimpse_postprocess:
	"""
	Glimpse does not keep track of the variant IDs in the INFO field.
	This rule adds them back in (assuming the variant representation is identical)
	"""
	input:
		glimpse = "{results}/external-calls/kage-glimpse/{sample}/{sample}_kage-glimpse_bi_genotyping_tmp_{chrom}.vcf.gz",
		vcf = "{results}/external-calls/panel-bi.vcf.gz"
	output:
		vcf = "{results}/external-calls/kage-glimpse/{sample}/{sample}_kage-glimpse_bi_genotyping_{chrom}.vcf.gz",
		tbi = "{results}/external-calls/kage-glimpse/{sample}/{sample}_kage-glimpse_bi_genotyping_{chrom}.vcf.gz.tbi"
	resources:
		mem_mb = 50000,
		walltime = "05:00:00"
	wildcard_constraints:
		chrom = "|".join(KAGE_CHROMOSOMES)
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		zcat {input.glimpse} | python3 workflow/scripts/add-ids.py {input.vcf} {wildcards.sample} | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""




rule external_kage_glimpse_concat_vcfs:
	"""
	Combine chromosome-specific KAGE and KAGE-Glimpse VCFs into a single VCF.
	"""
	input:
		vcfs = expand("{{results}}/external-calls/{{method}}/{{sample}}/{{sample}}_{{method}}_bi_genotyping_{chrom}.vcf.gz", chrom = KAGE_CHROMOSOMES)
	output:
		"{results}/external-calls/{method}/{sample}/{sample}_{method}_bi_genotyping.vcf.gz"
	log:
		"{results}/external-calls/{method}/{sample}/{sample}_{method}_bi_genotyping.log"
	benchmark:
		"{results}/external-calls/{method}/{sample}/{sample}_{method}_bi_genotyping.benchmark"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		method = "kage|kage-glimpse",
		chrom = "|".join(KAGE_CHROMOSOMES)
	threads: 15
	resources:
		mem_mb = 50000,
		walltime = "05:00:00"
	shell:
		"""
		bcftools concat -Ov --threads {threads} {input.vcfs} | python3 workflow/scripts/add_header_line.py | bgzip 1> {output} 2> {log}
		tabix -p vcf {output}
		"""
