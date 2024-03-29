#----- Merge Sequencing Tables -----#

rule seqtable_merge:
	"""
	Join sequence tables across all samples
	"""
	input:
		expand(os.path.join("results", "QC", "step_8", "clustered", "{sample}_seqtable.txt"), sample = SAMPLES)
	output:
		seqtable = os.path.join("results", "results", "seqtable.all"),
		tab2fa = os.path.join("results", "results", "seqtable.tab2fx")
	resources:
		mem_mb = 64000
	shell:
		"""
		{R}
		Rscript ./workflow/scripts/seqtable_merge.R
		"""
