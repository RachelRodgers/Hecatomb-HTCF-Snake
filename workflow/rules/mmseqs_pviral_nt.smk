#----- MMSeqs2 Query Probable Viral Unclassified Sequences (pviral_aa_unclassified_seqs.fasta) Against Refseq Virus NT UniVec Masked -----#

"""
Query pviral_aa_unclassified_seqs.fasta against
refseq_virus_nt_UniVec_masked/nt.fnaDB.
(All refseq viral gene sequences + nearest neighbors
with exact duplicates removed.)
"""

rule nt_create_querydb:
	input:
		os.path.join("results", "results", "mmseqs_aa_out", "pviral_aa_unclassified_seqs.fasta")
	output:
		os.path.join("results", "results", "mmseqs_nt_out", "seqtable_queryDB")
	shell:
		"""
		module load {MMSEQS}
		mmseqs createdb {input} {output} --dbtype 2
		"""

rule nt_search:
	"""
	Query the sequences from pviral aa unclassified against 
	refseq_virus_nt_UniVec_masked. This is a nucleotide search
	and will output an alignment.
	"""
	input:
		queryDB = os.path.join("results", "results", "mmseqs_nt_out", "seqtable_queryDB"),
		targetDB = NTTARGET
	params:
		alnDB = os.path.join("results", "results", "mmseqs_nt_out", "resultDB")
	output:
		idx = os.path.join("results", "results", "mmseqs_nt_out", "resultDB.index"),
		tmp = directory(os.path.join("results", "results", "mmseqs_nt_out", "tmp_nt"))
	threads: 16
	resources:
		mem_mb = 64000
	shell:
		"""
		module load {MMSEQS}
		mmseqs search {input.queryDB} {input.targetDB} \
			{params.alnDB} {output.tmp} \
			-a \
			-e 0.000001 \
			--search-type 3 \
			--cov-mode 2 \
			-c 0.95 \
			--threads {threads}
		"""

rule nt_extract_top_hit_from_search:
	input:
		idx = os.path.join("results", "results", "mmseqs_nt_out", "resultDB.index")
	params:
		resultDB = os.path.join("results", "results", "mmseqs_nt_out", "resultDB")
	output:
		bestResultDB = os.path.join("results", "results", "mmseqs_nt_out", "resultDB.firsthit")
	shell:
		"""
		module load {MMSEQS}
		mmseqs filterdb {params.resultDB} {output.bestResultDB} \
			--extract-lines 1 \
			--threads {threads}
		"""

rule nt_convert_top_hit_from_search:
	"""
	Convert best hit from search (resultDB.firsthit) to human readable
	"""
	input:
		queryDB = os.path.join("results", "results", "mmseqs_nt_out", "seqtable_queryDB"),
		targetDB = NTTARGET,
		alnDB = os.path.join("results", "results", "mmseqs_nt_out", "resultDB.firsthit")
	output:
		os.path.join("results", "results", "mmseqs_nt_out", "resultDB.firsthit.m8")
	resources:
		cpus = 1
		mem_mb = 1000
	shell:
		"""
		module load {MMSEQS}
		mmseqs convertalis \
			{input.queryDB} {input.targetDB} {input.alnDB} {output} \
			--threads {threads}
		"""
