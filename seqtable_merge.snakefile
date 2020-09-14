#----- Merge Sequencing Tables -----#

rule seqtable_merge:
	"""
	Join sequence tables across all samples
	"""
	input:
		expand(os.path.join("QC", "clustered", "{sample}_seqtable.txt"), sample = SAMPLES)
	output:
		seqtable = os.path.join("results", "seqtable.all"),
		tab2fa = os.path.join("results", "seqtable.tab2fx")
	shell:
		"""
		module load {R}
		./scripts/seqtable_merge.R
		"""