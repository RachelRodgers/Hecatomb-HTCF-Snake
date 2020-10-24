#----- Virome Assembly -----#

# Take QC'd reads from step 7 and perform metagenomic assembly with megahit.

rule assembly_digital_normalization:
	input:
		r1 = os.path.join("results", "QC", "step_7", PATTERN_R1 + ".s7.out.fastq"),
		r2 = os.path.join("results", "QC", "step_7", PATTERN_R2 + ".s7.out.fastq"),
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
			t={threads}
		"""

rule assembly_megahit:
	input:
		r1 = os.path.join("results", "assembly", PATTERN_R1 + ".norm.out.fastq"),
		r2 = os.path.join("results", "assembly", PATTERN_R2 + ".norm.out.fastq")
	params:
		prefix = "{sample}.mh"#,
		#outdir = directory(os.path.join("results", "assembly", "{sample}_megahit_out"))
	output:
		#fa = os.path.join("results", "assembly", "{sample}_megahit_out", "{sample}.mh.contigs.fa")
		outdir = directory(os.path.join("results", "assembly", "{sample}_megahit_out"))
	threads: 4
	shell:
		"""
		module load {MEGAHIT}
		megahit \
			-1 {input.r1} \
			-2 {input.r2} \
			-o {output.outdir} \
			--out-prefix {params.prefix}
		"""
