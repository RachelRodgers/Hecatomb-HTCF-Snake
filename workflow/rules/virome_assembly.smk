#----- Virome Assembly -----#

# Take QC'd reads from step 7 and perform metagenomic assembly with megahit.

rule assembly_kmer_stats:
	input:
		r1 = os.path.join("results", "QC", "step_7", PATTERN_R1 + ".s7.out.fastq"),
		r2 = os.path.join("results", "QC", "step_7", PATTERN_R2 + ".s7.out.fastq")
	output:
		os.path.join("results", "assembly", "stats", "{sample}_uniq_kmer_stats.txt")
	shell:
		"""
		module load {BBTOOLS}
		bbcountunique.sh \
			in={input.r1} \
			in2={input.r2} \
			interval=2500 \
			ow=t \
			{XMX} \
			out={output}
		"""

rule assembly_digital_normalization:
	input:
		r1 = os.path.join("results", "QC", "step_7", PATTERN_R1 + ".s7.out.fastq"),
		r2 = os.path.join("results", "QC", "step_7", PATTERN_R2 + ".s7.out.fastq")
		singletons = os.path.join("results", "QC", "step_7", "{sample}_singletons.s7.out.fastq")
	output:
		r1 = os.path.join("results", "assembly", PATTERN_R1 + ".norm.out.fastq"),
		r2 = os.path.join("results", "assembly", PATTERN_R2 + ".norm.out.fastq"),
		tossed = os.path.join("results", "assembly", "{sample}_tossed.norm.fastq"),
		hist = os.path.join("results", "assembly", "{sample}_norm.hist")
	threads: 4
	shell:
		"""
		module load {BBTOOLS}
		bbnorm.sh \
			in={input.r1} \
			in2={input.r2} \
			extra={input.singletons} \
			out={output.r1} \
			out2={output.r2} \
			outt={output.tossed} \
			hist={output.hist}
			target=20 mindepth=2 ow=t \
			{XMX} \
			t={threads}
		"""

rule assembly_megahit:
	input:
		r1 = os.path.join("results", "assembly", PATTERN_R1 + ".norm.out.fastq"),
		r2 = os.path.join("results", "assembly", PATTERN_R2 + ".norm.out.fastq")
	params:
		prefix = "{sample}_.mh"
	output:
		outdir = directory(os.path.join("results", "assembly", "{sample}_megahit_out"))
		fa = os.path.join("results", "assembly", "{sample}.mh.contigs.fa")
	threads: 4
	shell:
		"""
		module load {MEGAHIT}
		megahit \
			-1 {input.r1} \
			-2 {input.r2} \
			-o {output.outdir} \
			--out-prefix {params.prefix} \
			t {threads}
		"""

rule assembly_contig_copy:
	input:
		os.path.join("results", "assembly", "{sample}_megahit_out", "{sample}.mh.contigs.fa")
	output:
		os.path.join("results", "assembly", "megahit_contigs", "{sample}.mh.contigs.fa")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"cp {input} {output}"

rule assembly_quant_by_mapping:
	input:
		ref = os.path.join("results", "assembly", "megahit_contigs", "{sample}.mh.contigs.fa"),
		r1 = os.path.join("results", "QC", "step_7", PATTERN_R1 + ".s7.out.fastq"),
		r2 = os.path.join("results", "QC", "step_7", PATTERN_R2 + ".s7.out.fastq")
	output:
		os.path.join("results", "assembly", "megahit_contigs", "{sample}.aln.sam.gz")
	threads: 4
	shell:
		"""
		module load {BBTOOLS}
		bbmap.sh \
			ref={input.ref} \
			in={input.r1} \
			in2={input.r2} \
			out={output} \
			kfilter=22 subfilter=15 maxindel=80 ow=t \
			t={threads} \
			{XMX}
		"""

rule assembly_coverage:
	input:
		os.path.join("results", "assembly", "megahit_contigs", "{sample}.aln.sam.gz")
	output:
		os.path.join("results", "assembly", "megahit_contigs", "{sample}_cov.txt")
	shell:
		"""
		module load {BBTOOLS}
		pileup.sh in={input} out={output}
		"""

rule assembly_output_mapped_reads:
	input:
		os.path.join("results", "assembly", "megahit_contigs", "{sample}.aln.sam.gz")
	output:
		os.path.join("results", "assembly", "megahit_contigs", "{sample}_mapped.fastq")
	shell:
		"""
		module load {BBTOOLS}
		reformat.sh in={input} out={output} mappedonly
		"""

rule assembly_output_unmapped_reads:
	input:
		os.path.join("results", "assembly", "megahit_contigs", "{sample}.aln.sam.gz")
	output:
		os.path.join("results", "assembly", "megahit_contigs", "{sample}_mapped.fastq")
	shell:
		"""
		module load {BBTOOLS}
		reformat.sh in={input} out={output} unmappedonly
		"""