

rule histogram_extract_haplotype:
	"""
	Extract sequence from AGC file.
	"""
	input:
		AGC
	output:
		temp("{results}/assembly/{sample}.fa")
	conda:
		"..envs/histogram.yml"
	log:
		"{results}/{sample}/assembly/{sample}.log"
	resources:
		mem_mb = 40000
	shell:
		"""
		agc getset {input} {wildcards.sample} | bgzip > {output}
		samtools faidx {output}
		"""


rule histogram_run_biser:
	"""
	Detect SegDups with BISER.
	"""
	input:
		assembly = "{results}/{sample}/assembly/{sample}.fa",
		bed = lambda wildcards: BED[wildcards.sample],
	output:
		"{results}/{sample}/biser/{sample}_biser_segdups.bed"
	singularity:
		"workflow/container/biser.sif"
	log:
		"{results}/{sample}/biser/{sample}_segdups.log"
	resources:
		mem_mb = 100000
	threads: 5
	params:
		outname = "{results}/{sample}/biser/{sample}"
	shell:
		"""
		mkdir -p {wildcards.results}/{wildcards.sample}/biser/
		./workflow/scripts/biser.sh {input.bed} {input.assembly} {threads} {params.outname}
		"""


rule histogram_substract_repeats_from_segdups:
	"""
	Substract all repeat annotations overlapping segdups from
	repeat annotation BED to make sure that each region of the
	genome is either a repeat, a segdup or none of both.
	Then, create BED with repeats and segdups.
	"""
	input:
		repeatmasker_bed = lambda wildcards: BED[wildcards.sample],
		biser_bed =  "{results}/{sample}/biser/{sample}_biser_segdups.bed"
	output:
		repeats = "{results}/{sample}/plots/{sample}_repeats.bed",
		combined = "{results}/{sample}/plots/{sample}_repeats-segdups.bed"
	conda:
		"../envs/histogram.yml"
	shell:
		"""
		bedtools subtract -a {input.repeatmasker_bed} -b {input.biser_bed} | cut -f1,2,3,7 > {output.repeats}
		cat {output.repeats} {input.biser_bed} | bedtools sort > {output.combined}
		"""


rule histogram_plot_annotations:
	"""
	Plot histogram annotated by repeats and segdups.
	"""
	input:
		bed = "{results}/{sample}/plots/{sample}_repeats-segdups.bed",
		qvs = lambda wildcards: QVS[wildcards.sample]
	output:
		"{results}/{sample}/plots/{sample}_histogram.pdf"
	conda:
		"../envs/histogram.yml"
	params:
		outname = "{results}/{sample}/plots/{sample}"
	log:
		"{results}/{sample}/plots/{sample}_histogram.log"
	shell:
		"""
		mkdir -p {wildcards.results}/{wildcards.sample}/plots/
		./workflow/scripts/plot_qv_hist.sh {input.bed} {input.qvs} {params.outname} &> {log}
		"""
