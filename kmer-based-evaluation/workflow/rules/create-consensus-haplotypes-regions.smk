

rule subset_vcf:
	"""
	Extract only QV samples from VCF to reduce file size.
	This helps speeding up next steps.
	"""
	input:
		lambda wildcards: PHASED_VCFS[wildcards.callset]["vcf"]
	output:
		temp("{results}/haplotypes/{callset}/samples.vcf.gz")
	log:
		"{results}/haplotypes/{callset}/samples.log"
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
		bcftools view --threads {threads} -s {params.samples} {input} -O z -o {output} &> {log}
		tabix -p vcf {output}
		"""


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
	shell:
		"""
		bedtools maskfasta -fi {input.reference} -fo {output.fasta} -bed {input.bed} &> {log}
		samtools faidx {output.fasta}
		bgzip -c {output.fasta} > {output.fasta_gz}
		"""



rule regions_consensus_compute_consensus:
	"""
	Insert all variants into the reference genome to produce haplotypes.
	Remove all masked sequence regions.
	"""
	input:
		vcf = "{results}/haplotypes/{callset}/samples.vcf.gz",
		reference = "{results}/haplotypes/{callset}/{region}/genome/genome_{region}_masked.fasta",
		bed = "{results}/haplotypes/{callset}/{region}/bed/{region}_uncallable.bed"
	output:
		fasta_gz = temp("{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.fasta"),
		fasta_gz = temp("{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.fasta.gz"),
		temp_vcf = temp("{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.vcf.gz") 
	log:
		"{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.log"
	benchmark:
		"{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.benchmark.txt"
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		haplotype = "1|2"
	resources:
		mem_mb = 20000,
		walltime = "04:00:00"
	shell:
		"""
		
		bcftools consensus --sample {wildcards.sample} --haplotype {wildcards.haplotype} -e 'ALT~\"<.*>\"' -f {input.reference} {input.vcf} 2> {log} | python3 workflow/scripts/split_fasta.py -o {output.fasta} &> {log}
		bgzip -c {output.fasta} > {output.fasta_gz}
		"""


rule regions_consensus_compress_haplotypes:
	input:
		genomes = lambda wildcards: expand("{{results}}/haplotypes/{{callset}}/{{region}}_{{callset}}_{{region}}_{sample}_hap{haplotype}.fasta.gz", sample = SAMPLES[wildcards.callset], haplotype = ["1", "2"]),
		reference = lambda wildcards: PHASED_VCFS[wildcards.callset]["reference"]
	output:
		"{results}/haplotypes/{callset}_{region}.agc"
	log:
		"{results}/haplotypes/{callset}_{region}.log"
	wildcard_constraints:
		callset = "|".join([c for c in UNPHASED_VCFS.keys()] + [c for c in PHASED_VCFS.keys()])
	benchmark:
		"{results}/haplotypes/{callset}_{region}.benchmark.txt"
	conda:
		"../envs/agc.yaml"
	resources:
		mem_mb = 30000,
		walltime = "20:00:00"
	threads: 1 # 32
	shell:
		"""
		agc create {input.reference} {input.genomes} -o {output} -t {threads} &> {log}
		"""
