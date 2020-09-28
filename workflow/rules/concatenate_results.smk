#----- Concatenate Results -----#

"""
Adjust results to prepare for processing in R.
"""

rule fix_eukaryotic_virus_AA_annotation:
	input:
		os.path.join("results", "mmseqs_aa_checked_out", "viruses_checked_aa_table.tsv")
	output:
		os.path.join("results", "mmseqs_aa_checked_out", "viruses_checked_aa_table_edited.tsv")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"sed 's/uc_//g' {input} > {output}"

rule fix_eukaryotic_virus_NT_annotation:
	input:
		os.path.join("results", "mmseqs_nt_checked_out", "mmseqs_pviral_nt_checked_lineage.tsv")
	output:
		os.path.join("results", "mmseqs_nt_checked_out", "mmseqs_pviral_nt_checked_lineage_edited.tsv")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"tail -n +2 {input} | sed 's/NA/unknown/g' > {output}"

rule join_aa_and_nt_files:
	input:
		aa = os.path.join("results", "mmseqs_aa_checked_out", "viruses_checked_aa_table_edited.tsv"),
		nt = os.path.join("results", "mmseqs_nt_checked_out", "mmseqs_pviral_nt_checked_lineage_edited.tsv")
	output:
		os.path.join("results", "viruses_tax_table_tmp.tsv")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"cat {input.aa} {input.nt} | sort -n -k 1 > {output}"

rule fix_viruses_tax_table:
	input:
		os.path.join("results", "viruses_tax_table_tmp.tsv")
	output:
		os.path.join("results", "viruses_tax_table.tsv")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"sed '1iid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' {input} > {output}"

rule fix_phage_results_1:
	input:
		os.path.join("results", "mmseqs_aa_out", "phage_table.tsv")
	output:
		os.path.join("results", "phage_tax_table_tmp.tsv")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"sed 's/uc_//g' {input} > {output}"

rule fix_phage_results_2:
	input:
		os.path.join("results", "phage_tax_table_tmp.tsv")
	output:
		os.path.join("results", "phage_tax_table.tsv")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"sed '1iid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' {input} > {output}"

rule adjust_aa_aln:
	input:
		os.path.join("results", "mmseqs_aa_checked_out", "taxonomyResult.firsthit.m8")
	output:
		os.path.join("results", "mmseqs_aa_checked_out", "taxonomyResult.firsthit_edited.m8")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"sed '1iquery\ttarget\tpercent_id\talignment_length\tnum_mismatches\tnumber_gaps\tstart_query\tend_query\tstart_target\tend_target\te_value\tbit_score' {input} > {output}"

rule adjust_nt_aln:
	input:
		os.path.join("results", "mmseqs_nt_checked_out", "resultDB.firsthit.m8")
	output:
		os.path.join("results", "mmseqs_nt_checked_out", "resultDB.firsthit_edited.m8")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"sed '1iquery\ttarget\tpercent_id\talignment_length\tnum_mismatches\tnumber_gaps\tstart_query\tend_query\tstart_target\tend_target\te_value\tbit_score' {input} > {output}"

rule copy_aa_aln_to_results:
	input:
		os.path.join("results", "mmseqs_aa_checked_out", "taxonomyResult.firsthit_edited.m8")
	output:
		os.path.join("results", "aa.aln.m8")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"cp {input} {output}"

rule copy_nt_aln_to_results:
	input:
		os.path.join("results", "mmseqs_nt_checked_out", "resultDB.firsthit_edited.m8")
	output:
		os.path.join("results", "nt.aln.m8")
	resources:
		cpus=1,
		mem_mb=1000
	shell:
		"cp {input} {output}"
