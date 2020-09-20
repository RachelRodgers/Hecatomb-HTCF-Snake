#----- MMSeqs2 Query Probable Non-Phage Viral Sequences Against UniClust30 + Virus UniProt ProteinDB (Remove False Positives) -----#

"""
Query non-phage viral lineages pviral_virus_nt_seqs.fasta
(extracted from  mmseqs_pviral_nt_lineage.tsv) against
bac_virus_msaked/nt.fnaDB (UniClust 30 + Virus UniProt).
"""

rule ntcheck_extract_phage_lineages_from_pviralNT_tail:
	input:
		viruses = os.path.join("results", "mmseqs_nt_checked_out", "mmseqs_pviral_nt_lineage.tsv"),
		phage = PHAGE
	output:
		os.path.join("results", "mmseqs_nt_checked_out", "phage_nt_table.tsv")
	resources:
		cpus = 1, 
		mem_mb = 1000
	shell:
		"""
		tail -n+2 {input.viruses} | grep -f {input.phage} | sort -n -k1 > {output}
		"""

rule ntcheck_extract_phage_lineages_from_pviralNT_cut:
	input:
		os.path.join("results", "mmseqs_nt_checked_out", "phage_nt_table.tsv")
	output:
		os.path.join("results", "mmseqs_nt_checked_out", "phage_nt_table.list")
	resources:
		cpus = 1, 
		mem_mb = 1000
	shell:
		" cut -f1 {input} > {output}"

rule ntcheck_extract_phage_lineages_from_pviralNT_pullseq:
	input:
		seqtable = os.path.join("results", "seqtable.fasta"),
		list = os.path.join("results", "mmseqs_nt_checked_out", "phage_nt_table.list")
	output:
		os.path.join("results", "mmseqs_nt_checked_out", "phage_nt_seqs.fasta")
	shell:
		"""
		module load {PULLSEQ}
		pullseq -i {input.seqtable} -n {input.list} -l 5000 > {output}
		"""

rule ntcheck_extract_nonphage_viral_lineages_from_pviralNT_tail:
	input:
		viruses = os.path.join("results", "mmseqs_nt_checked_out", "mmseqs_pviral_nt_lineage.tsv"),
		phage = PHAGE
	output:
		os.path.join("results", "mmseqs_nt_checked_out", "pviral_virus_nt_table.tsv")
	resources:
		cpus = 1, 
		mem_mb = 1000
	shell:
		"""
		tail -n+2 {input.viruses} | grep -v -f {input.phage} | sort -n -k1 > {output}
		"""

rule ntcheck_extract_nonphage_viral_lineages_from_pviralNT_cut:
	input:
		os.path.join("results", "mmseqs_nt_checked_out", "pviral_virus_nt_table.tsv")
	output:
		os.path.join("results", "mmseqs_nt_checked_out", "pviral_virus_nt_table.list")
	resources:
		cpus = 1, 
		mem_mb = 1000
	shell:
		" cut -f1 {input} > {output}"

rule ntcheck_extract_nonphage_viral_lineages_from_pviralNT_pullseq:
	input:
		seqtable = os.path.join("results", "seqtable.fasta"),
		list = os.path.join("results", "mmseqs_nt_checked_out", "pviral_virus_nt_table.list")
	output:
		os.path.join("results", "mmseqs_nt_checked_out", "pviral_virus_nt_seqs.fasta")
	shell:
		"""
		module load {PULLSEQ}
		pullseq -i {input.seqtable} -n {input.list} -l 5000 > {output}
		"""

rule ntcheck_create_queryDB_from_nonphage_viral_lineages:
	"""
	Create a nucleotide query DB from the nonphage viral lineages
	"""
	input:
		os.path.join("results", "mmseqs_nt_checked_out", "pviral_virus_nt_seqs.fasta")
	output:
		os.path.join("results", "mmseqs_nt_checked_out", "seqtable_queryDB")
	shell:
		"""
		module load {MMSEQS}
		mmseqs createdb {input} {output} --dbtype 2
		"""

rule ntcheck_mmseqs_search:
	"""
	Query the sequences from pviral virus nt seqs against 
	bac_virus_masked. This is a nucleotide search
	and will output an alignment.
	"""
	input:
		queryDB = os.path.join("results", "mmseqs_nt_checked_out", "seqtable_queryDB"),
		targetDB = NTTARGETCHECK
	params:
		alnDB = os.path.join("results", "mmseqs_nt_checked_out", "resultDB")
	output:
		idx = os.path.join("results", "mmseqs_nt_checked_out", "resultDB.index"),
		tmp = directory(os.path.join("results", "mmseqs_nt_checked_out", "tmp_nt_check"))
	threads: 16
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

rule ntcheck_extract_top_hit_from_search:
	input:
		idx = os.path.join("results", "mmseqs_nt_checked_out", "resultDB.index")
	params:
		resultDB = os.path.join("results", "mmseqs_nt_checked_out", "resultDB")
	output:
		bestResultDB = os.path.join("results", "mmseqs_nt_checked_out", "resultDB.firsthit")
	shell:
		"""
		module load {MMSEQS}
		mmseqs filterdb {params.resultDB} {output.bestResultDB} \
			--extract-lines 1 \
			--threads {threads}
		"""

rule ntcheck_convert_top_hit_from_search:
	"""
	Convert best hit from search (resultDB.firsthit) to human readable
	"""
	input:
		queryDB = os.path.join("results", "mmseqs_nt_checked_out", "seqtable_queryDB"),
		targetDB = NTTARGETCHECK,
		alnDB = os.path.join("results", "mmseqs_nt_checked_out", "resultDB.firsthit")
	output:
		os.path.join("results", "mmseqs_nt_checked_out", "resultDB.firsthit.m8")
	shell:
		"""
		module load {MMSEQS}
		mmseqs convertalis \
			{input.queryDB} {input.targetDB} {input.alnDB} {output} \
			--threads {threads}
		"""
