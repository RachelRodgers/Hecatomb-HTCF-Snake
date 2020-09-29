#----- Annotate Probable Viral Non-Phage Sequences Search (Alignment) Results -----#

rule ntcheck_annotate:
	input:
		aln = os.path.join("results", "results", "mmseqs_nt_checked_out", "resultDB.firsthit.m8")
	output:
		os.path.join("results", "results", "mmseqs_nt_checked_out", "mmseqs_pviral_nt_checked_lineage.tsv")
	shell:
		"""
		module load {R}
		Rscript ./workflow/scripts/mmseqs_pviral_nt_check_annotate.R
		"""
