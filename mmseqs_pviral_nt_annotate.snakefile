#----- Annotate AA Unclassified Seqs NT Search (Alignment) Results -----#

rule mmseqs_pviral_nt_annotate:
	input:
		os.path.join("results", "mmseqs_nt_out", "resultDB.firsthit.m8")
	output:
		os.path.join("results", "mmseqs_nt_checked_out", "mmseqs_pviral_lineage.tsv")
	shell:
		"""
		modulel load {R}
		Rscript ./scripts/mmseqs_pviral_nt_annotate.R
		"""