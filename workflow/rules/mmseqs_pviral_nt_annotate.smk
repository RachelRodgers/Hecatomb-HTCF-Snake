#----- Annotate Probable Viral Unclassified Sequences Search (Alignment) Results -----#

rule nt_annotate:
	input:
		aln = os.path.join("results", "results", "mmseqs_nt_out", "resultDB.firsthit.m8")
	output:
		os.path.join("results", "results", "mmseqs_nt_checked_out", "mmseqs_pviral_nt_lineage.tsv")
	shell:
		"""
		{R}
		Rscript ./workflow/scripts/mmseqs_pviral_nt_annotate.R
		"""
