#----- Annotate AA Unclassified Seqs NT Search (Alignment) Results -----#

rule nt_annotate:
	input:
		aln = os.path.join("results", "mmseqs_nt_out", "resultDB.firsthit.m8")
	output:
		os.path.join("results", "mmseqs_nt_checked_out", "mmseqs_pviral_nt_lineage.tsv")
	shell:
		"""
		module load {R}
		Rscript ./scripts/mmseqs_pviral_nt_annotate.R
		"""
