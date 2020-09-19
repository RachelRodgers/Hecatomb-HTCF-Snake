#----- Annotate AA Unclassified Seqs NT Search (Alignment) Results -----#
rule nt_annotate_get_acc_path:
	run:
		# Write NCBI Accession path (for taxonomizr) to file so R can read in:
		ncbiAccPath = open("./taxonomizr_ncbi_accession_path.txt", "w")
		ncbiAccPath.write(NCBIACC)
		ncbiAccPath.close()

rule nt_annotate:
	input:
		aln = os.path.join("results", "mmseqs_nt_out", "resultDB.firsthit.m8"),
		path = os.path.join("taxonomizer_ncbi_accession_path.txt")
	output:
		os.path.join("results", "mmseqs_nt_checked_out", "mmseqs_pviral_nt_lineage.tsv")
	shell:
		"""
		module load {R}
		Rscript ./scripts/mmseqs_pviral_nt_annotate.R
		"""
