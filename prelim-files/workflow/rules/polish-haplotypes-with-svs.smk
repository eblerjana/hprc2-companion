
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
	conda:
		"../envs/whatshap.yml"
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
		cons = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa",
		cons_fai = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa.fai",
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


rule polish_extract_region_bam:
	"""
	Extract a region from the BAM file.
	"""
	input:
		chrom = "{results}/polishing/{callset}/contig-names/{sample}_{haplotype}/{chrom}.txt",
		bam = "{results}/polishing/{callset}/{tech}/{tech}_{sample}_{haplotype}.bam"
	output:
		"{results}/polishing/{callset}/{tech}/{tech}_{sample}_{haplotype}_{chrom}.bam"
	resources:
		mem_mb = 50000,
		walltime = "02:00:00"
	wildcard_constraints:
		tech = "illumina|ont",
		haplotype = "hap1|hap2"
	conda:
		"../envs/minimap.yml"
	log:
		"{results}/polishing/{callset}/{tech}/{tech}_{sample}_{haplotype}_{chrom}.log"
	benchmark:
		"{results}/polishing/{callset}/{tech}/{tech}_{sample}_{haplotype}_{chrom}.benchmark.txt"
	threads:
		10
	shell:
		"""
		samtools view -b --threads {threads} {input.bam} {wildcards.chrom} > {output}
		samtools index {output}
		"""


rule polish_call_small_variants:
	"""
	Call SNPs and indels from Illumina alignments using DeepVariant.
	"""
	input:
		bam = "{results}/polishing/{callset}/illumina/illumina_{sample}_{haplotype}_{chrom}.bam",
		ref = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa"
	output:
		vcf = "{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}_{chrom}.vcf.gz",
		gvcf = "{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}_{chrom}.g.vcf.gz"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	singularity:
		"workflow/container/deepvariant-1.8.0.sif"
	log:
		"{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}_{chrom}.log"
	benchmark:
		"{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}_{chrom}.benchmark.txt"
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
		 --regions {wildcards.chrom} \
		 --sample_name={wildcards.sample}-{wildcards.haplotype} \
		 --haploid_contigs="chrX,chrY,X,Y" &> {log}
		"""


rule polish_filter_calls:
	input:
		"{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}_{chrom}.vcf.gz"
	output:
		"{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}_{chrom}_filtered.vcf.gz"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	benchmark:
		"{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}_{chrom}_filtered.log"
	conda:
		"../envs/whatshap.yml"
	shell:
		"""
		bcftools view -f 'PASS' {input} | bgzip > {output}
		tabix -p vcf {output}
		"""


rule polish_genotype_variants:
	"""
	Re-genotype variants using WhatsHap and ONT reads.
	Only keep variants with matching genotypes between input and output.
	"""
	input: 
		bam = "{results}/polishing/{callset}/ont/ont_{sample}_{haplotype}_{chrom}.bam",
		vcf = "{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}_{chrom}_filtered.vcf.gz",
		reference = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa"
	output:
		vcf_filtered = "{results}/polishing/{callset}/whatshap-genotype/whatshap-genotype_{sample}_{haplotype}_{chrom}_filtered.vcf.gz",
		vcf_gz = "{results}/polishing/{callset}/whatshap-genotype/whatshap-genotype_{sample}_{haplotype}_{chrom}.vcf.gz", ## TODO make temp
		vcf = temp("{results}/polishing/{callset}/whatshap-genotype/whatshap-genotype_{sample}_{haplotype}_{chrom}.vcf")
	conda:
		"../envs/whatshap.yml"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	resources:
		mem_mb = 50000,
		walltime = "04:00:00"
	threads: 1
	log:
		"{results}/polishing/{callset}/whatshap-genotype/whatshap-genotype_{sample}_{haplotype}_{chrom}.log"
	benchmark:
		"{results}/polishing/{callset}/whatshap-genotype/whatshap-genotype_{sample}_{haplotype}_{chrom}.benchmark.txt"
	shell:
		" whatshap genotype --chromosome {wildcards.chrom} -o {output.vcf} --reference={input.reference} {input.vcf} {input.bam} &> {log} "
		" && "
		" cat {output.vcf} | python3 workflow/scripts/intersect-vcfs.py {input.vcf} | bgzip > {output.vcf_filtered} "
		" && "
		" tabix -p vcf {output.vcf_filtered} "
		" && "
		" bgzip -c {output.vcf} > {output.vcf_gz} "


rule polish_phase_variants:
	"""
	Phase variants using WhatsHap and ONT reads.
	"""
	input: 
		chrom = "{results}/polishing/{callset}/contig-names/{sample}_{haplotype}/{chrom}.txt",
		bam = "{results}/polishing/{callset}/ont/ont_{sample}_{haplotype}_{chrom}.bam",
		vcf = "{results}/polishing/{callset}/whatshap-genotype/whatshap-genotype_{sample}_{haplotype}_{chrom}_filtered.vcf.gz",
		reference = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa"
	output:
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
		" whatshap phase --chromosome {wildcards.chrom} -o {output.vcf} --reference={input.reference} {input.vcf} {input.bam} &> {log} "
		" && "
		" bgzip -c {output.vcf} > {output.vcf_gz} "
		" && "
		" tabix -p vcf {output.vcf_gz} "


rule polish_synchronize_variants:
	"""
	Synchronize phased blocks and filter VCF for haplotagging.
	"""
	input:
		"{results}/polishing/{callset}/whatshap/whatshap_{sample}_{haplotype}_{chrom}.vcf.gz"
	output:
		"{results}/polishing/{callset}/haplotag/haplotag_{sample}_{haplotype}_{chrom}.vcf.gz"
	conda:
		"../envs/whatshap.yml"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	benchmark:
		"{results}/polishing/{callset}/haplotag/haplotag_{sample}_{haplotype}_{chrom}.benchmark.txt"
	log:
		"{results}/polishing/{callset}/haplotag/haplotag_{sample}_{haplotype}_{chrom}.log"
	shell:
		"""
		python3 workflow/scripts/synchronize-phased-blocks.py {input} 2> {log} | bgzip > {output}
		tabix -p vcf {output}
		"""


rule polish_extract_haplotype_reads:
	"""
	Tag reads by haplotype and select haplotype of interest.
	"""
	input:
		vcf = "{results}/polishing/{callset}/haplotag/haplotag_{sample}_{haplotype}_{chrom}.vcf.gz",
		bam = "{results}/polishing/{callset}/ont/ont_{sample}_{haplotype}_{chrom}.bam",
		reference = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa"
	output:
		tagged = temp("{results}/polishing/{callset}/haplotag/haplotag_{sample}_{haplotype}_{chrom}.bam"),
		split = "{results}/polishing/{callset}/haplotag/split_{sample}_{haplotype}_{chrom}.bam",
		hap1 = temp("{results}/polishing/{callset}/haplotag/split_{sample}_{haplotype}_{chrom}_tmp1.bam"),
		hap2 = temp("{results}/polishing/{callset}/haplotag/split_{sample}_{haplotype}_{chrom}_tmp2.bam"),
		tsv = "{results}/polishing/{callset}/haplotag/haplotag_{sample}_{haplotype}_{chrom}.tsv"
	conda:
		"../envs/whatshap.yml"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	log:
		haplotag = "{results}/polishing/{callset}/haplotag/haplotag_{sample}_{haplotype}_{chrom}_haplotag.log",
		split = "{results}/polishing/{callset}/haplotag/haplotag_{sample}_{haplotype}_{chrom}_split.log"
	benchmark:
		"{results}/polishing/{callset}/haplotag/haplotag_{sample}_{haplotype}_{chrom}.benchmark.txt"
	shell:	
		" whatshap haplotag -o {output.tagged} --reference {input.reference} {input.vcf} {input.bam} --output-haplotag-list {output.tsv} &> {log.haplotag} "
		" && "
		" whatshap split --output-h1 {output.hap1} --output-h2 {output.hap2} {input.bam} {output.tsv} &> {log.split} "
		" && " 
		" samtools sort {output.hap1} | samtools markdup - {output.split} "
		" && " 
		" samtools index {output.split} "


rule polish_call_svs:
	"""
	Call SVs using Sniffles2.
	Keep only homozygous calls and remove BND calls.
	"""
	input:
		bam = "{results}/polishing/{callset}/haplotag/split_{sample}_{haplotype}_{chrom}.bam",
		reference = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa"
	output:
		vcf = "{results}/polishing/{callset}/sniffles/sniffles_{sample}_{haplotype}_{chrom}.vcf.gz",
		filtered = "{results}/polishing/{callset}/sniffles/sniffles_{sample}_{haplotype}_{chrom}_filtered.vcf.gz"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	resources:
		mem_mb = 50000,
		walltime = "02:00:00"
	conda:
		"../envs/sniffles.yml"
	benchmark:
		"{results}/polishing/{callset}/sniffles/sniffles_{sample}_{haplotype}_{chrom}.benchmark.txt"
	log:
		"{results}/polishing/{callset}/sniffles/sniffles_{sample}_{haplotype}_{chrom}.log"
	shell:
		" sniffles -i {input.bam} -v {output.vcf} --reference {input.reference} --sample-id {wildcards.sample}-{wildcards.haplotype}  &> {log} "
		" && "
		" bcftools view  -f 'PASS' --min-ac 2 {output.vcf} | bcftools filter -i \'FMT/DR<1 & FMT/DV>5\' | grep -v \"BND\" | bgzip > {output.filtered} "
		" && " 
		" tabix -p vcf {output.filtered} "


rule polish_detect_errors:
	"""
	Detect errors in the consensus haplotypes.
	"""
	input:
		small = "{results}/polishing/{callset}/haplotag/haplotag_{sample}_{haplotype}_{chrom}.vcf.gz",
		svs = "{results}/polishing/{callset}/sniffles/sniffles_{sample}_{haplotype}_{chrom}_filtered.vcf.gz"
	output:
		vcf = "{results}/polishing/{callset}/errors/errors_{sample}_{haplotype}_{chrom}.vcf.gz",
		vcf_small = temp("{results}/polishing/{callset}/haplotag/haplotag_{sample}_{haplotype}_{chrom}_tmp.vcf.gz"),
		bed = temp("{results}/polishing/{callset}/sniffles/sniffles_{sample}_{haplotype}_{chrom}.bed")
	wildcard_constraints:
		haplotype = "hap1|hap2"
	conda:
		"../envs/whatshap.yml"
	log:
		"{results}/polishing/{callset}/errors/errors_{sample}_{haplotype}_{chrom}.log"
	benchmark:
		"{results}/polishing/{callset}/errors/errors_{sample}_{haplotype}_{chrom}.benchmark.txt"
	shell:
		" bcftools query -f \'%CHROM\\t%POS0\\t%END\\n\' {input.svs}  > {output.bed} "
		" && "
		" bedtools subtract -header -a {input.small} -b {output.bed} | python3 workflow/scripts/find_errors.py -window-width {WINDOW_WIDTH} 2> {log} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' | bgzip > {output.vcf_small} "
		" && "
		" tabix -p vcf {output.vcf_small} "
		" && "
		" bcftools concat -a {output.vcf_small} {input.svs} -Oz -o {output.vcf} "
		" && "
		" tabix -p vcf {output.vcf} "



def aggregate_concat_phased_variants(wildcards):
	checkpoint_output = checkpoints.polish_determine_contig_names.get(**wildcards).output[0]
	return expand("{results}/polishing/{callset}/errors/errors_{sample}_{haplotype}_{chrom}.vcf.gz",
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
		"{results}/polishing/{callset}/errors/errors_{sample}_{haplotype}.vcf.gz"
	conda:
		"../envs/whatshap.yml"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	resources:
		mem_mb = 50000,
		walltime = "02:00:00"
	threads: 15
	log:
		"{results}/polishing/{callset}/errors/errors_whatshap_{sample}_{haplotype}.log"
	benchmark:
		"{results}/polishing/{callset}/errors/errors_whatshap_{sample}_{haplotype}.benchmark.txt"
	shell:
		"""
		bcftools concat -o {output} -Oz --threads {threads} {input.vcfs} &> {log}
		bcftools index {output}
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
		bcftools consensus --haplotype A -f {input.haplotype}  -e \'ALT~\"<.*>\"\' {input.errors} 2> {log}  > {output}
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


