#----- Contaminant Removal -----#

# Remove non-biological sequences, low-qual bases, host and (100% ID) bacterial sequences from virome-sequencied libraries.

rule clumpify:
	"""
	Step 0: Clumpify & deduplicate reads
	"""
	input:
		r1 = os.path.join(READDIR + "/renamed/", PATTERN_R1 + ".fastq.gz"),
		r2 = os.path.join(READDIR + "/renamed/", PATTERN_R2 + ".fastq.gz")
	output:
		r1 = os.path.join("results", "clumped", PATTERN_R1 + ".clumped.fastq.gz"),
		r2 = os.path.join("results", "clumped", PATTERN_R2 + ".clumped.fastq.gz")
	threads: 8
	shell:
		"""
		module load {BBTOOLS}
		clumpify.sh \
			in={input.r1} \
			in2={input.r2} \
			out={output.r1} \
			out2={output.r2} \
			reorder=a \
			ow=t \
			t={threads}
		"""

rule remove_leftmost_primerB:
	"""
	Step 1: Remove leftmost primerB
	"""
	input:
		r1 = os.path.join("results", "clumped", PATTERN_R1 + ".clumped.fastq.gz"),
		r2 = os.path.join("results", "clumped", PATTERN_R2 + ".clumped.fastq.gz"),
		primers = os.path.join(CONPATH, "primerB.fa")
	output:
		r1 = temp(os.path.join("results", "QC", "step_1", PATTERN_R1 + ".s1.out.fastq")),
		r2 = temp(os.path.join("results", "QC", "step_1", PATTERN_R2 + ".s1.out.fastq")),
                stats = os.path.join("results", "QC", "step_1", "{sample}.s1.stats")
	threads: 8
	shell:
		"""
		module load {BBTOOLS}
		bbduk.sh \
			in={input.r1} \
			in2={input.r2} \
			ref={input.primers} \
			out={output.r1} \
			out2={output.r2} \
			stats={output.stats} \
			k=16 hdist=1 mink=11 ktrim=l restrictleft=20 \
			removeifeitherbad=f \
			trimpolya=10 \
			ordered=t \
			rcomp=f \
			ow=t \
			t={threads}
		"""

rule remove_3prime_contaminant:
	"""
	Step 2: Remove 3' read-thru contaminant"
	"""
	input:
		r1 = os.path.join("results", "QC", "step_1", PATTERN_R1 + ".s1.out.fastq"),
		r2 = os.path.join("results", "QC", "step_1", PATTERN_R2 + ".s1.out.fastq"),
		primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
	output:
		r1 = temp(os.path.join("results", "QC", "step_2", PATTERN_R1 + ".s2.out.fastq")),
		r2 = temp(os.path.join("results", "QC", "step_2", PATTERN_R2 + ".s2.out.fastq")),
		stats = os.path.join("results", "QC", "step_2", "{sample}.s2.stats")
	threads: 8
	shell:
		"""
		module load {BBTOOLS}
		bbduk.sh \
			in={input.r1} \
			in2={input.r2} \
			ref={input.primers} \
			out={output.r1} \
			out2={output.r2} \
			stats={output.stats} \
			k=16 hdist=1 mink=11 ktrim=r \
			removeifeitherbad=f \
			ordered=t \
			rcomp=f \
			ow=t \
			threads={threads}
		"""

rule remove_primer_free_adapter:
	"""
	Step 3: Remove primer free adapter
	"""
	input:
		r1 = os.path.join("results", "QC", "step_2", PATTERN_R1 + ".s2.out.fastq"),
		r2 = os.path.join("results", "QC", "step_2", PATTERN_R2 + ".s2.out.fastq"),
		primers = os.path.join(CONPATH, "nebnext_adapters.fa")
	output:
		r1 = temp(os.path.join("results", "QC", "step_3", PATTERN_R1 + ".s3.out.fastq")),
		r2 = temp(os.path.join("results", "QC", "step_3", PATTERN_R2 + ".s3.out.fastq")),
		stats = os.path.join("results", "QC", "step_3", "{sample}.s3.stats")
	threads: 8
	shell:
		"""
		module load {BBTOOLS}
		bbduk.sh \
			in={input.r1} \
			in2={input.r2} \
			ref={input.primers} \
			out={output.r1} \
			out2={output.r2} \
			stats={output.stats} \
			k=16 hdist=1 mink=10 ktrim=r \
			removeifeitherbad=f \
			ordered=t \
			rcomp=t \
			ow=t \
			threads={threads}
		"""

rule remove_adapter_free_primer:
	"""
	Step 4: Remove adapter-free primer
	"""
	input:
		r1 = os.path.join("results", "QC", "step_3", PATTERN_R1 + ".s3.out.fastq"),
		r2 = os.path.join("results", "QC", "step_3", PATTERN_R2 + ".s3.out.fastq"),
		primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
	output:
		r1 = temp(os.path.join("results", "QC", "step_4", PATTERN_R1 + ".s4.out.fastq")),
		r2 = temp(os.path.join("results", "QC", "step_4", PATTERN_R2 + ".s4.out.fastq")),
		stats = os.path.join("results", "QC", "step_4", "{sample}.s4.stats")
	threads: 8
	shell:
		"""
		module load {BBTOOLS}
		bbduk.sh \
			in={input.r1} \
			in2={input.r2} \
			ref={input.primers} \
			out={output.r1} \
			out2={output.r2} \
			stats={output.stats} \
			k=16 hdist=0 \
			removeifeitherbad=f \
			ordered=t \
			rcomp=t \
			ow=t \
			t={threads}
		"""

rule remove_vector_contamination:
	"""
	Step 5: Vector contamination removal (PhiX + NCBI UniVecDB)
	"""
	input:
		r1 = os.path.join("results", "QC", "step_4", PATTERN_R1 + ".s4.out.fastq"),
		r2 = os.path.join("results", "QC", "step_4", PATTERN_R2 + ".s4.out.fastq"),
		primers = os.path.join(CONPATH, "vector_contaminats.fa.gz")
	output:
		r1 = temp(os.path.join("results", "QC", "step_5", PATTERN_R1 + ".s5.out.fastq")),
		r2 = temp(os.path.join("results", "QC", "step_5", PATTERN_R2 + ".s5.out.fastq")),
		stats = os.path.join("results", "QC", "step_5", "{sample}.s5.stats")
	threads: 8
	shell:
		"""
		module load {BBTOOLS}
		bbduk.sh \
			in={input.r1} \
			in2={input.r2} \
			ref={input.primers} \
			out={output.r1} \
			out2={output.r2} \
			stats={output.stats} \
			k=31 hammingdistance=1 \
			ordered=t \
			ow=t \
			t={threads}
		"""

rule host_removal:
	"""
	Step 6a: Host removal
	"""
	input:
		r1 = os.path.join("results", "QC", "step_5", PATTERN_R1 + ".s5.out.fastq"),
		r2 = os.path.join("results", "QC", "step_5", PATTERN_R2 + ".s5.out.fastq"),
		reference = HOSTPATH
	output:
		unmapped = temp(os.path.join("results", "QC", "step_6", "{sample}_unmapped.s6.out.fastq")),
		mapped = temp(os.path.join("results", "QC", "step_6", "{sample}_hostmapped.s6.out.fastq"))
	threads: 8
	resources:
		mem_mb=50000
	shell:
		"""
		module load {BBTOOLS}
		bbmap.sh \
			in={input.r1} \
			in2={input.r2} \
			outu={output.unmapped} \
			outm={output.mapped} \
			semiperfectmode=t \
			quickmatch fast \
			ordered=t \
			path={input.reference} \
			{XMX} \
			ow=t \
			t={threads}
		"""

rule repair:
	"""
	Step 6b: Repair the paired ends
	"""
	input:
		unmapped=os.path.join("results", "QC", "step_6", "{sample}_unmapped.s6.out.fastq")
	output:
		r1 = temp(os.path.join("results", "QC", "step_6", PATTERN_R1 + ".s6.out.fastq")),
		r2 = temp(os.path.join("results", "QC", "step_6", PATTERN_R2 + ".s6.out.fastq"))
	shell:
		"""
		module load {BBTOOLS}
		repair.sh \
			in={input.unmapped} \
			out={output.r1} \
			out2={output.r2} \
			ow=t
		"""

rule trim_low_quality:
	"""
	Step 7a: Trim low quality bases
	"""
	input:
		r1 = os.path.join("results", "QC", "step_6", PATTERN_R1 + ".s6.out.fastq"),
		r2 = os.path.join("results", "QC", "step_6", PATTERN_R2 + ".s6.out.fastq")
	output:
		r1 = os.path.join("results", "QC", "step_7", PATTERN_R1 + ".s7.out.fastq"),
		r2 = os.path.join("results", "QC", "step_7", PATTERN_R2 + ".s7.out.fastq"),
		singletons = os.path.join("results", "QC", "step_7", "{sample}_singletons.s7.out.fastq"),
		stats = os.path.join("results", "QC", "step_7", "{sample}.s7.stats")
	threads: 8
	shell:
		"""
		module load {BBTOOLS}
		bbduk.sh \
			in={input.r1} \
			in2={input.r2} \
			out={output.r1} \
			out2={output.r2} \
			outs={output.singletons} \
			stats={output.stats} \
			qtrim=r trimq=20 \
			maxns=2 minlength=50 \
			ordered=t \
			ow=t \
			threads={threads}
		"""
rule get_r1_singletons:
	"""
	Step 7b: Split R1 singletons
	"""
	input:
		singletons = os.path.join("results", "QC", "step_7", "{sample}_singletons.s7.out.fastq")
	output:
		r1singletons = os.path.join("results", "QC", "step_7", "{sample}_singletons_R1.out.fastq")
	shell:
		"""
		grep -A 3 '1:N:' {input.singletons} | sed '/^--$/d' > {output.r1singletons}
		"""

rule get_r2_singletons:
	"""
	Step 7c: Split R2 singletons
	"""
	input:
		singletons = os.path.join("results", "QC", "step_7", "{sample}.singletons.s7.out.fastq")
	output:
		r2singletons = os.path.join("results", "QC", "step_7", "{sample}_singletons_R2.out.fastq")
	shell:
		"""
		grep -A 3 '2:N:' {input.singletons} | sed '/^--$/d' > {output.r2singletons}
		"""

rule concat_r1:
	"""
	Step 7d: Concatenate R1 reads & singletons
	"""
	input:
		r1 = os.path.join("results", "QC", "step_7", PATTERN_R1 + ".s7.out.fastq"),
		r1singletons = os.path.join("results", "QC", "step_7", "{sample}_singletons_R1.out.fastq")
	output:
		r1combo = os.path.join("results", "QC", "step_7", PATTERN_R1 + ".s7.combined.out.fastq")
	shell:
		"""
		cat {input.r1} {input.r1singletons} > {output.r1combo} 
		"""

rule concat_r2:
	"""
	Step 7e: Concatenate R2 reads & singletons
	"""
	input:
                r2 = os.path.join("results", "QC", "step_7", PATTERN_R2 + ".s7.out.fastq"),
                r2singletons = os.path.join("results", "QC", "step_7", "{sample}_singletons_R2.out.fastq")
	output:
		r2combo = os.path.join("results", "QC", "step_7", PATTERN_R2 + ".s7.combined.out.fastq")
	shell:
		"""
		cat {input.r2} {input.r2singletons} > {output.r2combo}
		"""

rule remove_bacteria:
	"""
	Step 8: Remove bacterial contaminants reserving viral and ambiguous sequences
	"""
	input:
		r1combo = os.path.join("results", "QC", "step_7", PATTERN_R1 + ".s7.combined.out.fastq"),
		reference = BACPATH 
	output:
		mapped = os.path.join("results", "QC", "step_8", "{sample}_bacterial.fastq"),
		unmapped = os.path.join("results", "QC", "step_8", "{sample}_viral_amb.fastq"),
		scafstats = os.path.join("results", "QC", "step_8", "{sample}_scafstats.txt")
	threads: 16
	resources:
		mem_mb=50000
	shell:
		"""
		module load {BBTOOLS}
		bbmap.sh \
			in={input.r1combo} \
			path={input.reference} \
			outm={output.mapped} \
			outu={output.unmapped} \
			semiperfectmode=t \
			quickmatch fast \
			scafstats={output.scafstats} \
			ordered=t \
			{XMX} \
			ow=t \
			t={threads}
		"""
