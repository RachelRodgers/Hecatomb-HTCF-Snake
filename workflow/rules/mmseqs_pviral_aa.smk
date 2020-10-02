#----- MMSeqs2 Query Sequencing Table Against UniProt Viral Proteins Clustered at 99% (Virus Uniprot)-----#
# Extract phage (phage_tax_table.tsv), non-phage (viruses_seqs.fasta), and unclassified (pviral_aa_unclassified_seqs.fasta)

"""
Query sample sequences against Virus Uniprot DB.
Generate a phage sequence table & alignment,
non-phage virus sequence table & alignment,
and an unclassified table.
"""

rule aa_convert_seqtable_to_fasta:
	"""
	Convert seqtable.tab2fx to fasta so a sequence db can be generated from it
	"""
	input:
		os.path.join("results", "results", "seqtable.tab2fx")
	output:
		os.path.join("results", "results", "seqtable.fasta")
	shell:
		"""
		module load {SEQKIT}
		seqkit tab2fx {input} -o {output} -w 5000 -o {output}
		"""

rule aa_create_querydb_from_seqtable:
	"""
	Create a sequence db from the seqtable made from all samples
	"""
	input:
		os.path.join("results", "results", "seqtable.fasta")
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "seqtable_queryDB")
	shell:
		"""
		module load {MMSEQS}
		mmseqs createdb {input} {output} --dont-shuffle 0 --dbtype 0
		"""

rule aa_taxonomy_search_alignment:
	"""
	Query the sequences in seqtable_queryDB (from the seqtable generated from all samples)
	against the virus uniprotDB to generate taxonomic assignments to the sequences using LCA.
	This is a translated search and will output an alignment.
	"""
	input:
		queryDB = os.path.join("results", "results", "mmseqs_aa_out", "seqtable_queryDB"),
		targetDB = AATARGET
	params:
		taxaDB = os.path.join("results", "results", "mmseqs_aa_out", "tax_search_alignment", "taxonomyResult")
	output:
		tmp = directory(os.path.join("results", "results", "mmseqs_aa_out", "tax_search_alignment", "tmp_aa")),
		idx = os.path.join("results", "results", "mmseqs_aa_out", "tax_search_alignment", "taxonomyResult.index")
	threads: 16
	resources:
		mem_mb = 64000
	shell:
		"""
		module load {MMSEQS}
		mmseqs taxonomy \
			{input.queryDB} {input.targetDB} {params.taxaDB} {output.tmp} \
			-a \
			--start-sens 1 \
			--sens-steps 3 \
			-s 7 \
			--search-type 2 \
			--tax-output-mode 1 \
			--threads {threads}
		"""

rule aa_convert_taxonomy_alignment_results_to_m8:
	"""
	Convert the alignment results DB (taxonomyResult) to a human-readable format
	"""
	input:
		queryDB = os.path.join("results", "results", "mmseqs_aa_out", "seqtable_queryDB"),
		targetDB = AATARGET,
		idx = os.path.join("results", "results", "mmseqs_aa_out", "tax_search_alignment", "taxonomyResult.index")
	params:
		alnDB = os.path.join("results", "results", "mmseqs_aa_out", "tax_search_alignment", "taxonomyResult")
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "aln.m8")
	resources:
		mem_mb = 4000
	shell:
		"""
		module load {MMSEQS}
		mmseqs convertalis \
			{input.queryDB} {input.targetDB} {params.alnDB} {output} \
			--format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln"
		"""

rule aa_taxonomy_search_lca:
	"""
	Query the sequences in seqtable_queryDB (from the seqtable generated from all samples)
	against the virus uniprotDB to generate taxonomic assignments to the sequences using LCA.
	This is a translated search and will output LCA (taxonomy).
	"""
	input:
		queryDB = os.path.join("results", "results", "mmseqs_aa_out", "seqtable_queryDB"),
		targetDB = AATARGET
	params:
		taxaDB = os.path.join("results", "results", "mmseqs_aa_out", "lcaDB")
	output:
		tmp = directory(os.path.join("results", "results", "mmseqs_aa_out", "tmp_aa")),
		idx = os.path.join("results", "results", "mmseqs_aa_out", "lcaDB.index")
	threads: 16
	resources:
		mem_mb = 64000
	shell:
		"""
		module load {MMSEQS}
		mmseqs taxonomy \
			{input.queryDB} {input.targetDB} {params.taxaDB} {output.tmp} \
			-a \
			--start-sens 1 \
			--sens-steps 3 \
			-s 7 \
			--search-type 2 \
            --tax-lineage true \
			--lca-ranks "superkingdom:phylum:class:order:family:genus:species"
		"""

rule aa_convert_taxonomy_lca_results_to_tsv:
	"""
	Create a TSV formatted taxonomy table from the taxonomy LCA output
	"""
	input:
		queryDB = os.path.join("results", "results", "mmseqs_aa_out", "seqtable_queryDB"),
		idx = os.path.join("results", "results", "mmseqs_aa_out", "lcaDB.index")
	params:
		resultDB = os.path.join("results", "results", "mmseqs_aa_out", "lcaDB")
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "taxonomyResult.tsv")
	threads: 8
	resources:
		mem_mb = 4000
	shell:
		"""
		module load {MMSEQS}
		mmseqs createtsv {input.queryDB} {params.resultDB} {output} \
			--threads {threads}
		"""

rule aa_extract_all_potential_viruses:
	"""
	Extract all potential viral sequences from taxonomyResult.tsv
	"""
	input:
		os.path.join("results", "results", "mmseqs_aa_out", "taxonomyResult.tsv")
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "all_viruses_table.tsv")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"""
		grep 'Viruses:' {input} | cut -f1,5 | sed 's/phi14:2/phi14_2/g' | sed 's/:/\t/g' | sort -n -k1 > {output}
		"""

rule aa_extract_phage_lineages_for_R_grep:
	"""
	Extract phage lineages from all_viruses_table.tsv and generate
	taxonomy table for import into R as a phyloseq object
	"""
	input:
		viruses = os.path.join("results", "results", "mmseqs_aa_out", "all_viruses_table.tsv"),
		phagetax = PHAGE
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "phage_table.tsv")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"grep -f {input.phagetax} {input.viruses} > {output}"

rule aa_extract_phage_lineages_for_R_cut:
	input:
		os.path.join("results", "results", "mmseqs_aa_out", "phage_table.tsv")
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "phage_seqs.list")
	resources:
            cpus=1,
            mem_mb=1000
	shell:
		"cut -f1 {input} > {output}"

rule aa_extract_phage_lineages_for_R_pullseq:
	input:
		seqtable = os.path.join("results", "results", "seqtable.fasta"),
		list = os.path.join("results", "results", "mmseqs_aa_out", "phage_seqs.list")
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "phage_seqs.fasta")
	shell:
		"""
		module load {PULLSEQ} 
        pullseq -i {input.seqtable} -n {input.list} -l 5000 > {output}
		"""

rule aa_extract_phage_lineages_for_R_seqkit:
	input:
		os.path.join("results", "results", "mmseqs_aa_out", "phage_seqs.fasta")
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "phage_seqs.fx2tab")
	shell:
		"""
		module load {SEQKIT}
        seqkit fx2tab {input} > {output}
		"""

rule aa_extract_phage_lineages_for_R_join:
	input:
		seqs = os.path.join("results", "results", "mmseqs_aa_out", "phage_seqs.fx2tab"),
		table = os.path.join("results", "results", "mmseqs_aa_out", "phage_table.tsv")
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "phage_tax_table.tsv")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"""
		join {input.seqs} {input.table} | \
			awk -F ' ' '{{ print $2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",$7,"\t",$8,"\t",$9 }}' |
			sed '1isequence\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' > {output}
	"""

rule aa_extract_nonphage_viral_lineages_for_R_grep:
	input:
		viruses = os.path.join("results", "results", "mmseqs_aa_out", "all_viruses_table.tsv"),
		phagetax = PHAGE,
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "viruses_table.tsv")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"grep -v -f {input.phagetax} {input.viruses} > {output}"

rule aa_extract_nonphage_viral_lineages_for_R_cut:
	input:
		os.path.join("results", "results", "mmseqs_aa_out", "viruses_table.tsv")
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "viruses_seqs.list")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"cut -f1 {input} > {output}"

rule aa_extract_nonphage_viral_lineages_for_R_pullseq:
	input:
		seqtable = os.path.join("results", "results", "seqtable.fasta"),
		list = os.path.join("results", "results", "mmseqs_aa_out", "viruses_seqs.list")
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "viruses_seqs.fasta")
	shell:
		"""
		module load {PULLSEQ}
		pullseq -i {input.seqtable} -n {input.list} -l 5000 > {output}
		"""

rule aa_extract_unclassified_grep:
	input:
		os.path.join("results", "results", "mmseqs_aa_out", "taxonomyResult.tsv")
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "pviral_aa_unclassified_seqs.list")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"""
		grep -v 'Viruses:' {input} | cut -f1,5 | sed 's/:/\t/g' | \
			sort -n -k1 > {output}
		"""

rule aa_extract_unclassified_pullseq:
	input:
		seqtable = os.path.join("results", "results", "seqtable.fasta"),
		list = os.path.join("results", "results", "mmseqs_aa_out", "pviral_aa_unclassified_seqs.list")
	output:
		os.path.join("results", "results", "mmseqs_aa_out", "pviral_aa_unclassified_seqs.fasta")
	shell:
		"""
		module load {PULLSEQ}
		pullseq -i {input.seqtable} -n {input.list} -l 5000 > {output}
		"""
