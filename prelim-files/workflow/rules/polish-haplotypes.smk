
rule polish_fasta_to_fastq:
	"""
	If fasta files are provided, add constant
	quality scores
	"""
	input:
		lambda wildcards: ILLUMINA[wildcards.sample]
	output:
		"{results}/polishing/fastqs/{sample}.fastq.gz"
	benchmark:
		"{results}/polishing/fastqs/{sample}.benchmark.txt"
	resources:
		walltime = "10:00:00"
	shell:
		"""
		zcat {input} | awk 'BEGIN {{RS = \">\" ; FS = \"\\n\"}} NR > 1 {{print \"@\"$1\"\\n\"$2\"\\n+\"$1\"\\n\"gensub(/./, \":\", \"g\", $2)}}' | bgzip > {output}
		"""


rule polish_extract_sequence:
	"""
	Extract consensus haplotype from agc archieve
	and index it.
	"""
	input:
		lambda wildcards: AGC[wildcards.callset]
	output:
		cons = temp("{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa"),
		cons_fai = temp("{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa.fai"),
		bwt = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa.bwt"
	resources:
		mem_mb = 30000
	conda:
		"../envs/bwa.yml"
	benchmark:
		"{results}/polishing/{callset}/consensus/{sample}_{haplotype}.benchmark.txt"
	params:
		name = lambda wildcards: CONSENSUS[wildcards.callset][(wildcards.sample, wildcards.haplotype)]
	shell:
		"""
		agc getset {input} {params.name} > {output.cons}
		bwa index {output.cons}
		samtools faidx {output.cons}
		"""


rule polish_align_illumina:
	"""
	Align Illumina reads to consensus haplotype.
	"""
	input:
		haplotype = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa",
		bwt = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa.bwt",
		reads = lambda wildcards: ILLUMINA[wildcards.sample] if ILLUMINA[wildcards.sample].endswith(('fastq.gz', 'fq.gz')) else "{results}/polishing/fastqs/{sample}.fastq.gz"
	output:
		"{results}/polishing/{callset}/illumina/illumina_{sample}_{haplotype}.bam"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	resources:
		mem_mb = 90000,
		walltime = "25:00:00" 
	conda:
		"../envs/bwa.yml"
	log:
		bwa = "{results}/polishing/{callset}/illumina/illumina_{sample}_{haplotype}.log",
		sam = "{results}/polishing/{callset}/illumina/illumina_{sample}_{haplotype}.sam.log"
	benchmark:
		"{results}/polishing/{callset}/illumina/illumina_{sample}_{haplotype}.benchmark.txt"
	threads: 24
	params:
		sam_threads = 8
	shell:
		" bwa mem -t {threads} -p -M {input.haplotype} {input.reads} 2> {log.bwa}"
		" | "
		" samtools view -bS --threads {params.sam_threads} " 
		" | "
		" samtools sort --threads {params.sam_threads} -T {wildcards.callset}_{wildcards.sample}_{wildcards.haplotype}_bwa -o {output} 2> {log.sam} "
		" && " 
		" samtools index --threads {threads} {output}"


rule polish_align_ont:
	"""
	Align ONT reads to consensus haplotype.
	"""
	input:
		haplotype = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa",
		reads = lambda wildcards: ONT[wildcards.sample],
		cram_ref = CRAM_REF
	output:
		"{results}/polishing/{callset}/ont/ont_{sample}_{haplotype}.bam"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	resources:
		mem_mb = 30000,
		walltime = "25:00:00"
	conda:
		"../envs/minimap.yml"
	log:
		mm2 = "{results}/polishing/{callset}/ont/ont_{sample}_{haplotype}.log",
		sam = "{results}/polishing/{callset}/ont/ont_{sample}_{haplotype}.sam.log"
	benchmark:
		 "{results}/polishing/{callset}/ont/ont_{sample}_{haplotype}.benchmark.txt"
	threads: 24
	params:
		readgroup = lambda wildcards:  (f'"@RG\\tID:{wildcards.sample}-{wildcards.haplotype}'f'\\tSM:{wildcards.sample}-{wildcards.haplotype}"'),
		sam_threads = 8
	shell:
		"samtools fastq --reference {input.cram_ref} {input.reads} | minimap2 -ax map-ont -Y --MD --eqx -t {threads} -R {params.readgroup} --secondary=no {input.haplotype} - 2> {log.mm2} "
		" | "
		" samtools view -bS --threads {params.sam_threads} "
		" | "
		" samtools sort --threads {params.sam_threads} -T {wildcards.sample}_{wildcards.haplotype}_mm2 -o {output} 2> {log.sam}"
		" && "
		" samtools index --threads {threads} {output}"


rule polish_call_small_variants:
	"""
	Call SNPs and indels from Illumina alignments using DeepVariant.
	"""
	input:
		bam = "{results}/polishing/{callset}/illumina/illumina_{sample}_{haplotype}.bam",
		ref = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa"
	output:
		vcf = "{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}.vcf.gz",
		gvcf = "{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}.g.vcf.gz"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	singularity:
		"workflow/container/deepvariant-1.8.0.sif"
	log:
		"{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}.log"
	benchmark:
		"{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}.benchmark.txt"
	threads: 24
	resources:
		mem_mb = 50000,
		walltime = "07:00:00"
	shell:
		"""
		 /opt/deepvariant/bin/run_deepvariant \
		 --model_type=WGS \
		 --ref={input.ref} \
		 --reads={input.bam} \
		 --output_vcf={output.vcf} \
		 --output_gvcf={output.gvcf} \
		 --num_shards={threads} \
		 --vcf_stats_report=true \
		 --disable_small_model=true \
		 --sample_name={wildcards.sample}-{wildcards.haplotype} \
		 --haploid_contigs="chrX,chrY,X,Y" &> {log}
		"""

rule polish_filter_calls:
	input:
		"{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}.vcf.gz"
	output:
		"{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}_filtered.vcf.gz"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	benchmark:
		"{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}_filtered.log"
	conda:
		"../envs/whatshap.yml"
	shell:
		"""
		bcftools view -f 'PASS' {input} | bcftools norm -m- -Oz -o {output}
		tabix -p vcf {output}
		"""



checkpoint polish_determine_contig_names:
	"""
	Determine contig names present in consensus haplotypes.
	This is necessary to run whatshap separately on each contig.
	"""
	input:
		"{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa.fai"
	output:
		temp(directory("{results}/polishing/{callset}/contig-names/{sample}_{haplotype}"))
	wildcard_constraints:
		haplotype = "hap1|hap2"
	params:
		outprefix = "{results}/polishing/{callset}/contig-names/{sample}_{haplotype}/"
	shell:
		"""
		mkdir -p {output}
		python3 workflow/scripts/determine-contig-names.py {input} {params.outprefix}
		"""


rule polish_phase_variants:
	"""
	Phase variants using WhatsHap and ONT reads.
	Extract
	"""
	input: 
		chrom = "{results}/polishing/{callset}/contig-names/{sample}_{haplotype}/{chrom}.txt",
		bam = "{results}/polishing/{callset}/ont/ont_{sample}_{haplotype}.bam",
		vcf = "{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}_filtered.vcf.gz",
		reference = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa"
	output:
		region_vcf = temp("{results}/polishing/{callset}/whatshap/input_{sample}_{haplotype}_{chrom}.vcf.gz"),
		region_tbi = temp("{results}/polishing/{callset}/whatshap/input_{sample}_{haplotype}_{chrom}.vcf.gz.tbi"),
		vcf_gz = "{results}/polishing/{callset}/whatshap/whatshap_{sample}_{haplotype}_{chrom}.vcf.gz",
		vcf = temp("{results}/polishing/{callset}/whatshap/whatshap_{sample}_{haplotype}_{chrom}.vcf")
	conda:
		"../envs/whatshap.yml"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	resources:
		mem_mb = 50000,
		walltime = "04:00:00"
	threads: 1
	log:
		"{results}/polishing/{callset}/whatshap/whatshap_{sample}_{haplotype}_{chrom}.log"
	benchmark:
		"{results}/polishing/{callset}/whatshap/whatshap_{sample}_{haplotype}_{chrom}.benchmark.txt"
	shell:
		"""
		bcftools view -r {wildcards.chrom} {input.vcf} -o {output.region_vcf} -Oz
		tabix -p vcf {output.region_vcf}
		whatshap phase --chromosome {wildcards.chrom} -o {output.vcf} --reference={input.reference} {output.region_vcf} {input.bam} &> {log}
		bgzip -c {output.vcf} > {output.vcf_gz}
		tabix -p vcf {output.vcf_gz}
		"""


def aggregate_concat_phased_variants(wildcards):
	checkpoint_output = checkpoints.polish_determine_contig_names.get(**wildcards).output[0]
	prefix = "whatshap_" + wildcards.sample + "_" + wildcards.haplotype
	return expand("{results}/polishing/{callset}/whatshap/whatshap_{sample}_{haplotype}_{chrom}.vcf.gz",
			results=wildcards.results,
			callset=wildcards.callset,
			sample=wildcards.sample,
			haplotype=wildcards.haplotype,
			chrom=glob_wildcards(os.path.join(checkpoint_output, "{chrom}.txt")).chrom)


rule polish_concat_phased_variants:
	"""
	Concat chromosome-wise phased VCFs into a single VCF.
	"""
	input:
		vcfs = aggregate_concat_phased_variants
	output:
		"{results}/polishing/{callset}/whatshap/phased_whatshap_{sample}_{haplotype}.vcf.gz"
	conda:
		"../envs/whatshap.yml"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	resources:
		mem_mb = 50000,
		walltime = "02:00:00"
	threads: 15
	log:
		"{results}/polishing/{callset}/whatshap/phased_whatshap_{sample}_{haplotype}.log"
	benchmark:
		"{results}/polishing/{callset}/whatshap/phased_whatshap_{sample}_{haplotype}.benchmark.txt"
	shell:
		"""
		bcftools concat -o {output} -Oz --threads {threads} {input.vcfs} &> {log}
		bcftools index {output}
		"""


rule polish_detect_errors:
	"""
	Detect errors in the consensus haplotypes.
	"""
	input:
		"{results}/polishing/{callset}/whatshap/phased_whatshap_{sample}_{haplotype}.vcf.gz"
	output:
		"{results}/polishing/{callset}/errors/errors_{sample}_{haplotype}.vcf.gz"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	log:
		"{results}/polishing/{callset}/errors/errors_{sample}_{haplotype}.log"
	benchmark:
		"{results}/polishing/{callset}/errors/errors_{sample}_{haplotype}.benchmark.txt"
	shell:
		"""
		zcat {input} | python3 workflow/scripts/find_errors.py -window-width {WINDOW_WIDTH} 2> {log} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' | bgzip > {output}
		tabix -p vcf {output}
		"""


rule polish_correct_errors:
	"""
	Correct errors of the consensus haplotypes.
	"""
	input:
		errors = "{results}/polishing/{callset}/errors/errors_{sample}_{haplotype}.vcf.gz",
		haplotype = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa"
	output:
		"{results}/polishing/{callset}/polished/{sample}_{haplotype}.fa"
	log:
		"{results}/polishing/{callset}/polished/{sample}_{haplotype}.log"
	benchmark:
		"{results}/polishing/{callset}/polished/{sample}_{haplotype}.benchmark.txt"
	conda:
		"../envs/whatshap.yml"
	shell:
		"""
		bcftools consensus --haplotype A -f {input.haplotype} {input.errors} 2> {log}  > {output}
		"""

rule polish_compress_haplotypes:
	input:
		expand("{{results}}/polishing/{{callset}}/polished/{sample}_{haplotype}.fa", sample = SAMPLES, haplotype = ["hap1", "hap2"])
	output:
		"{results}/polishing/{callset}/polished/{callset}_polished.agc"
	benchmark:
		"{results}/polishing/{callset}/polished/{callset}_polished.benchmark.txt"
	log:
		"{results}/polishing/{callset}/polished/{callset}_polished.log"
	conda:
		"../envs/bwa.yml"
	threads:
		24
	resources:
		mem_mb = 100000,
		walltime = "02:00:00"
	shell:
		"""
		agc create {input} -o {output} -t {threads} &> {log}
		"""



rule polish_count_variants:
	"""
	Compute variant statistics.
	"""
	input:
		"{filename}.vcf.gz"
	output:
		"{filename}.stats"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	conda:
		"../envs/rtg.yml"
	shell:
		"""
		rtg vcfstats {input} &> {output}
		"""


