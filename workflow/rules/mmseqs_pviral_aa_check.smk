#----- MMSeqs2 Query Probable Viral Non-Phage Sequences (viruses_seqs.fasta) Against UniClust30 ProteinDB (Remove False Positives; Uni Plus Virus)-----#

"""
Query probable viral sequences (viruses_seqs.fasta) against UniClust30 DB.
This step is to remove false positives ("check").
The end result of these steps is a checked viral sequence table and a checked
viral alignment table.
UniClust 30 is all UniProtKB entires clustered at 30%ID concatenated to
Virus UniProt entires clustered at 99%.
"""

rule aacheck_create_querydb_from_viruses_seqs:
	"""
	Create a sequence db from the viruses_seqs.fasta file (probably non-phage viral sequences)
	"""
	input:
		os.path.join("results", "results", "mmseqs_aa_out", "viruses_seqs.fasta")
	output:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "viral_seqs_queryDB")
	shell:
		"""
		{MMSEQS} createdb \
			{input} {output} --dbtype 0
		"""

rule aacheck_taxonomy_search_alignment:
	"""
	Query the sequences in viral_seqs_queryDB (generated from all probable non-phage
	viral sequences) against the UniClust 30 protein DB to generate taxonomic
	assignments to the sequences using LCA.  This is a translated search and will output 
	an alignment.
	"""
	input:
		queryDB = os.path.join("results", "results", "mmseqs_aa_checked_out", "viral_seqs_queryDB"),
		targetDB = AATARGETCHECK
	params:
		taxaDB = os.path.join("results", "results", "mmseqs_aa_checked_out", "tax_search_alignment", "taxonomyResult")
	output:
		tmp = directory(os.path.join("results", "results", "mmseqs_aa_checked_out", "tax_search_alignment", "tmp_aa_checked")),
		idx = os.path.join("results", "results", "mmseqs_aa_checked_out", "tax_search_alignment", "taxonomyResult.index")
	threads: 16
	resources:
		mem_mb = 64000
	shell:
		"""
		{MMSEQS} taxonomy \
			{input.queryDB} {input.targetDB} {params.taxaDB} {output.tmp} \
			-a \
			-s 7 \
			--search-type 2 \
			--tax-output-mode 1 \
			--threads {threads}
		"""

rule aacheck_convert_taxonomy_result_to_m8:
	"""
	Convert the alignment results DB (taxonomyResult) to a human-readable format
	"""
	input:
		queryDB = os.path.join("results", "results", "mmseqs_aa_checked_out", "viral_seqs_queryDB"),
		targetDB = AATARGETCHECK,
		idx = os.path.join("results", "results", "mmseqs_aa_checked_out", "tax_search_alignment", "taxonomyResult.index")
	params:
		alnDB = os.path.join("results", "results", "mmseqs_aa_checked_out", "tax_search_alignment", "taxonomyResult")
	output:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "aln.m8")
	resources:
		mem_mb = 4000
	shell:
		"""
		{MMSEQS} convertalis \
			{input.queryDB} {input.targetDB} {params.alnDB} {output} \
			--format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln" \
			--threads 16
		"""

rule aacheck_extract_best_hit_from_taxonomy_result:
	input:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "tax_search_alignment", "taxonomyResult.index")
	params:
		resultDB = os.path.join("results", "results", "mmseqs_aa_checked_out", "tax_search_alignment", "taxonomyResult"),
		bestResultDB = os.path.join("results", "results", "mmseqs_aa_checked_out", "taxonomyResult.firsthit")
	output:
		bestResultDB = os.path.join("results", "results", "mmseqs_aa_checked_out", "taxonomyResult.firsthit.dbtype")
	shell:
		"""
		{MMSEQS} filterdb \
			{params.resultDB} {params.bestResultDB} --extract-lines 1
		"""

rule aacheck_convert_best_hit:
	"""
	Convert best hit from taxonomyResult (taxonomyResult.firsthit) to human readable
	"""
	input:
		queryDB = os.path.join("results", "results", "mmseqs_aa_checked_out", "viral_seqs_queryDB"),
		targetDB = AATARGETCHECK,
		idx = os.path.join("results", "results", "mmseqs_aa_checked_out", "taxonomyResult.firsthit.dbtype")
	params:
		alignmentDB = os.path.join("results", "results", "mmseqs_aa_checked_out", "taxonomyResult.firsthit")
	output:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "taxonomyResult.firsthit.m8")
	threads: 4
	shell:
		"""
		{MMSEQS} convertalis \
			{input.queryDB} {input.targetDB} {params.alignmentDB} {output} --threads {threads}
		"""

rule aacheck_taxonomy_search_lca:
	"""
	Query the sequences in viral_seqs_queryDB (generated from all probable non-phage
	viral sequences) against the UniClust 30 protein DB to generate taxonomic
	assignments to the sequences using LCA.  This is a translated search will output 
	LCA (taxonomy).
	"""
	input:
		queryDB = os.path.join("results", "results", "mmseqs_aa_checked_out", "viral_seqs_queryDB"),
		targetDB = AATARGETCHECK
	params:
		taxaDB = os.path.join("results", "results", "mmseqs_aa_checked_out", "tax_search_lca", "lcaDB")
	output:
		tmp = directory(os.path.join("results", "results", "mmseqs_aa_checked_out", "tax_search_lca", "tmp_aa_checked")),
		idx = os.path.join("results", "results", "mmseqs_aa_checked_out", "tax_search_lca", "lcaDB.index")
	threads: 16
	resources:
		mem_mb = 64000
	shell:
		"""
		{MMSEQS} taxonomy \
			{input.queryDB} {input.targetDB} {params.taxaDB} {output.tmp} \
			-a \
			-s 7 \
			--search-type 2 \
			--tax-lineage true \
			--lca-ranks "superkingdom,phylum,class,order,family,genus,species" \
			--threads {threads}
		"""

rule aacheck_create_taxonomy_table_from_lca:
	"""
	Create a TSV formatted taxonomy table from the LCA output
	"""
	input:
		queryDB = os.path.join("results", "results", "mmseqs_aa_checked_out", "viral_seqs_queryDB"),
		idx = os.path.join("results", "results", "mmseqs_aa_checked_out", "tax_search_lca", "lcaDB.index")
	params:
		resultDB = os.path.join("results", "results", "mmseqs_aa_checked_out", "tax_search_lca", "lcaDB")
	output:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "taxonomyResult.tsv")
	threads: 4
	resources:
		mem_mb = 4000
	shell:
		"""
		{MMSEQS} createtsv \
			{input.queryDB} {params.resultDB} {output} \
			--threads {threads}
		"""

rule aacheck_extract_nonphage_viral_lineages_for_R_grep:
	"""
	Extract non-phage viral lineages from taxonomyResult.tsv and generate
	taxonomy table for import into R as a phyloseq object
	"""
	input:
		viruses = os.path.join("results", "results", "mmseqs_aa_checked_out", "taxonomyResult.tsv"),
		phagetax = PHAGE
	output:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "viruses_checked_aa_table.tsv")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"""
		grep -v 'Bacteria:' {input.viruses} | grep 'Viruses:' | grep -v -f {input.phagetax} | cut -f1,5 | sed 's/:/\t/g' | sort -n -k1 > {output}
		"""

rule aacheck_extract_nonphage_viral_lineages_for_R_cut:
	input:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "viruses_checked_aa_table.tsv")
	output:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "viruses_checked_aa_seqs.list")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"cut -f1 {input} > {output}"

rule aacheck_extract_nonphage_viral_lineages_for_R_pullseq:
	input:
		seqtable = os.path.join("results", "results", "seqtable.fasta"),
		list = os.path.join("results", "results", "mmseqs_aa_checked_out", "viruses_checked_aa_seqs.list")
	output:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "viruses_checked_aa_seqs.fasta")
	shell:
		"""
		{PULLSEQ} -i {input.seqtable} -n {input.list} -l 5000 > {output}
		"""

rule aacheck_extract_nonphage_viral_lineages_for_R_seqkit:
	input:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "viruses_checked_aa_seqs.fasta")
	output:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "viruses_checked_aa_seqs.fx2tab")
	shell:
		"""
		{SEQKIT}
		seqkit fx2tab {input} > {output}
		"""

rule aacheck_extract_nonphage_viral_lineages_for_R_join:
	input:
		fx2tab = os.path.join("results", "results", "mmseqs_aa_checked_out", "viruses_checked_aa_seqs.fx2tab"),
		tsv = os.path.join("results", "results", "mmseqs_aa_checked_out", "viruses_checked_aa_table.tsv")
	output:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "viruses_checked_aa_tax_table.tsv")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"""
		join {input.fx2tab} {input.tsv} | \
			awk -F ' ' '{{ print $2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",$7,"\t",$8,"\t",$9 }}' | \
			sed '1isequence\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' > \
			{output}
		"""

rule aacheck_extract_unclassified_lineages_grep:
	input:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "taxonomyResult.tsv")
	output:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "unclassified_checked_aa_seqs.list")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"grep -v 'Viruses:' {input} | cut -f1,5 | sed 's/:/\t/g' | sort -n -k1 > {output}"

rule aacheck_extract_unclassified_lineages_pullseq:
	input:
		seqtable = os.path.join("results", "results", "seqtable.fasta"),
		list = os.path.join("results", "results", "mmseqs_aa_checked_out", "unclassified_checked_aa_seqs.list")
	output:
		os.path.join("results", "results", "mmseqs_aa_checked_out", "unclassified_checked_aa_seqs.fasta")
	shell:
		"""
		{PULLSEQ} -i {input.seqtable} -n {input.list} -l 5000 > {output}
		"""
