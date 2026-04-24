

rule regions_compute_complement:
	"""
	This produces BED files defining all regions outside of the
	region of interest.
	"""
	input:
		bed = lambda wildcards: PHASED_VCFS[wildcards.callset]["regions"][wildcards.region],
		fai = lambda wildcards: PHASED_VCFS[wildcards.callset]["reference"] + ".fai"
	output:
		tmp_bed = temp("{results}/haplotypes/{callset}/{region}/bed/{region}_uncallable_tmp.bed"),
		tmp_fai = temp("{results}/haplotypes/{callset}/{region}/bed/{region}_genome.fai"),
		bed = "{results}/haplotypes/{callset}/{region}/bed/{region}_uncallable.bed"
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		region = "|".join(REGIONS)
	shell:
		"""
		cat {input.bed} | sort -k1,1 -k2,2n > {output.tmp_bed}
		sort -k1,1 -k2,2n {input.fai} > {output.tmp_fai} 
		bedtools complement -i {output.tmp_bed} -g {output.tmp_fai} > {output.bed}
		"""


rule regions_mask_genome:
	"""
	Set all bases outside of the desired regions to N.
	"""
	input:
		reference = lambda wildcards: PHASED_VCFS[wildcards.callset]["reference"],
		bed = "{results}/haplotypes/{callset}/{region}/bed/{region}_uncallable.bed"
	output:
		fasta = temp("{results}/haplotypes/{callset}/{region}/genome/genome_{region}_masked.fasta"),
		fasta_gz = temp("{results}/haplotypes/{callset}/{region}/genome/genome_{region}_masked.fasta.gz")
	log:
		"{results}/haplotypes/{callset}/{region}/genome/genome_{region}_masked.log"
	conda:
		"../envs/shapeit.yaml"
	resources:
		mem_total_mb = 10000
	wildcard_constraints:
		region = "|".join(REGIONS)
	shell:
		"""
		bedtools maskfasta -fi {input.reference} -fo {output.fasta} -bed {input.bed} &> {log}
		samtools faidx {output.fasta}
		bgzip -c {output.fasta} > {output.fasta_gz}
		"""


rule regions_subset_vcf:
	"""
	Extract only QV samples from VCF to reduce file size.
	This helps speeding up next steps.
	"""
	input:
		vcf = lambda wildcards: PHASED_VCFS[wildcards.callset]["vcf"],
		bed = "{results}/haplotypes/{callset}/{region}/bed/{region}_uncallable.bed"
	output:
		tmp = temp("{results}/haplotypes/{callset}/samples_{region}_tmp.vcf.gz"),
		vcf = temp("{results}/haplotypes/{callset}/samples_{region}.vcf.gz")
	log:
		"{results}/haplotypes/{callset}/samples_{region}.log"
	conda:
		"../envs/shapeit.yaml"
	threads:
		24
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 4
	params:
		samples = lambda wildcards: ",".join(SAMPLES[wildcards.callset])
	shell:
		"""
		bcftools view --threads {threads} -s {params.samples} {input.vcf} -O z -o {output.tmp} &> {log}
		tabix -p vcf {output.tmp}
		bedtools subtract -header -A -a {output.tmp} -b {input.bed} | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


rule regions_consensus_compute_consensus:
	"""
	Insert all variants into the reference genome to produce haplotypes.
	Remove all masked sequence regions.
	"""
	input:
		vcf = "{results}/haplotypes/{callset}/samples_{region}.vcf.gz",
		reference = "{results}/haplotypes/{callset}/{region}/genome/genome_{region}_masked.fasta"
	output:
		fasta = temp("{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.fasta"),
		fasta_gz = temp("{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.fasta.gz")
	log:
		"{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.log"
	benchmark:
		"{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.benchmark.txt"
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		haplotype = "1|2",
		region = "|".join(REGIONS)
	resources:
		mem_mb = 20000,
		walltime = "04:00:00"
	shell:
		"""	
		bcftools consensus --sample {wildcards.sample} --haplotype {wildcards.haplotype} -e 'ALT~\"<.*>\"' -f {input.reference} {input.vcf} 2> {log} | python3 workflow/scripts/split_fasta.py -o {output.fasta} &> {log}
		bgzip -c {output.fasta} > {output.fasta_gz}
		"""
