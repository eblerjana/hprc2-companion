

rule consensus_compute_consensus:
	"""
	Insert all variants into the reference genome to produce haplotypes.
	"""
	input:
		vcf = lambda wildcards: PHASED_VCFS[wildcards.callset]["vcf"],
		reference = lambda wildcards: PHASED_VCFS[wildcards.callset]["reference"]
	output:
		temp("{results}/haplotypes/{callset}/all/{callset}_all_{sample}_hap{haplotype}.fasta.gz")
	log:
		"{results}/haplotypes/{callset}/all/{callset}_all_{sample}_hap{haplotype}.log"
	benchmark:
		"{results}/haplotypes/{callset}/all/{callset}_all_{sample}_hap{haplotype}.benchmark.txt"
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		haplotype = "1|2"
	resources:
		mem_mb = 20000,
		walltime = "04:00:00"
	shell:
		"""
		bcftools consensus --sample {wildcards.sample} --haplotype {wildcards.haplotype} -e 'ALT~\"<.*>\"' -f {input.reference} {input.vcf} 2> {log} | bgzip > {output}
		"""


rule consensus_compress_haplotypes:
	input:
		genomes = lambda wildcards: expand("{{results}}/haplotypes/{{callset}}/{{region}}/{{callset}}_{{region}}_{sample}_hap{haplotype}.fasta.gz", sample = SAMPLES[wildcards.callset], haplotype = ["1", "2"]),
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
