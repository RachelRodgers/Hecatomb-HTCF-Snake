#----- Cluster Count -----#

# Dereplicate and count sequences

rule remove_exact_duplicates:
	"""
	Step 1: Remove exact duplicates
	"""
	input:
		os.path.join("results", "QC", "step_8", "{sample}_viral_amb.fastq")
	output:
		os.path.join("results", "QC", "step_8", "clustered", PATTERN_R1 + ".s8.deduped.out.fastq")
	threads: 4
	shell:
		"""
		bash {BBTOOLS} \
		dedupe.sh \
			in={input} \
			ow=t \
			out={output} \
			ac=f \
			{XMX} \
			t={threads}
		"""

rule dereplicate:
	"""
	Step 2a: Dereplicate
	"""
	input:
		os.path.join("results", "QC", "step_8", "clustered", PATTERN_R1 + ".s8.deduped.out.fastq")
	output:
		fa = os.path.join("results", "QC", "step_8", "clustered", "{sample}_best.fasta"),
		stats = os.path.join("results", "QC", "step_8", "clustered", "{sample}_stats.txt")
	threads: 8
	resources:
		mem_mb = 250000
	shell:
		"""
		bash {BBTOOLS} \
		dedupe.sh \
			in={input} \
			ow=t \
			s=4 \
			rnc=t \
			pbr=t \
			csf={output.stats} \
			out={output.fa} \
			{XMX} \
			t={threads}
		"""

rule reformat:
	"""
	Step 3a: Reformat
	"""
	input:
		os.path.join("results", "QC", "step_8", "clustered", "{sample}_best.fasta")
	output:
		out = os.path.join("results", "QC", "step_8", "clustered", "{sample}_reformatted.fasta")
	shell:
		"""
		bash {BBTOOLS} \
		reformat.sh \
			in={input} \
			out={output} \
			deleteinput=t \
			fastawrap=0 \
			ow=t \
			{XMX}
		""" 

rule extract_sequences:
	"""
	Step 3b: Extract sequence lines from sample_reformatted.fasta
	"""
	input:
		os.path.join("results", "QC", "step_8", "clustered", "{sample}_reformatted.fasta")
	output:
		os.path.join("results", "QC", "step_8", "clustered", "{sample}_seqs.txt")
	resources:
		mem_mb = 4000
	shell:
		"grep -v '>' {input} | sed '1i sequence' > {output}"

rule extract_counts:
	"""
	Step 3c: Extract counts from dedupe stats file
	"""
	input:
		os.path.join("results", "QC", "step_8", "clustered", "{sample}_stats.txt")
	output:
		os.path.join("results", "QC", "step_8", "clustered", "{sample}_counts.txt")
	resources:
		mem_mb = 4000
	shell:
		"cut -f 2 {input} | sed '1s/size/{wildcards.sample}/' > {output}"

rule create_sequence_table:
	"""
	Step 3d: Create sequence table
	"""
	input:
		seqs = os.path.join("results", "QC", "step_8", "clustered", "{sample}_seqs.txt"),
		counts = os.path.join("results", "QC", "step_8", "clustered", "{sample}_counts.txt")
	output:
		os.path.join("results", "QC", "step_8", "clustered", "{sample}_seqtable.txt")
	resources:
		mem_mb = 4000
	shell:
		"paste {input.seqs} {input.counts} > {output}"
