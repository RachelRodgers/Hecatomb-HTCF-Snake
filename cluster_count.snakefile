#----- Cluster Count -----#

# Dereplicate and count sequences

rule remove_exact_duplicates:
	"""
	Step 9: Remove exact duplicates
	"""
	input:
		os.path.join("QC", "step_8", "{sample}_viral_amb.fastq")
	output:
		os.path.join("QC", "step_9", PATTERN_R1 + ".s9.deduped.out.fastq")
	threads: 8
	resources: 
		mem_mb = 50000
	shell:
		"""
		module load {BBTOOLS}
		dedupe.sh \
			in={input} \
			out={output} \
			ow=t ac=f \
			{XMX} \
			t={threads}
		"""

rule dereplicate:
	"""
	Step 10: Dereplicate
	"""
	input:
		os.path.join("QC", "step_9", PATTERN_R1 + ".s9.deduped.out.fastq")
	output:
		fa = os.path.join("QC", "step_10", "{sample}_best.fasta"),
		stats = os.path.join("QC", "step_10", "{sample}_stats.txt")
	threads: 8
	resources:
		mem_mb = 50000
	shell:
		"""
		module load {BBTOOLS}
		dedupe.sh \
			in={input} \
			out={output.fa} \
			csf={output.stats} \
			ow=t s=4 rnc=t pbr=t \
			{XMX} \
			t={threads}
		"""

rule extract_sequence_counts:
	"""
	Step 11: Extract sequences and counts for seqtable (count table)
	"""
	input:
		os.path.join("QC", "step_10", "{sample}_best.fasta")
	output:
		out = os.path.join("QC", "step_11", "{sample}_reformatted.fasta")
	resources:
		mem_mb = 50000

	shell:
		"""
		module load {BBTOOLS}
		reformat.sh \
			in={input} \
			out={output} \
			deleteinput=t fastawrap=0 ow=t \
			{XMX}
		""" 

rule extract_sequences:
	"""
	Step 12: Extract sequence lines from sample_reformatted.fasta
	"""
	input:
		os.path.join("QC", "step_11", "{sample}_reformatted.fasta")
	output:
		os.path.join("QC", "clustered", "{sample}_seqs.txt")
	shell:
		"grep -v '>' {input} | sed '1i sequence' > {output}"

rule extract_counts:
	"""
	Step 13: Extract counts from dedupe stats file
	"""
	input:
		os.path.join("QC", "step_10", "{sample}_stats.txt")
	output:
		os.path.join("QC", "clustered", "{sample}_counts.txt")
	shell:
		"cut -f 2 {input} | sed '1s/size/{wildcards.sample}/' > {output}"

rule create_sequence_table:
	"""
	Step 14: Create sequence table
	"""
	input:
		seqs = os.path.join("QC", "clustered", "{sample}_seqs.txt"),
		counts = os.path.join("QC", "clustered", "{sample}_counts.txt")
	output:
		os.path.join("QC", "clustered", "{sample}_seqtable.txt")
	shell:
		"paste {input.seqs} {input.counts} > {output}"