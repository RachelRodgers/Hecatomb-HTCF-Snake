"""
Master Snakefile for running the Hecatomb pipeline on HTCF.

Rachel Rodgers, Sep 2020
""" 

import os
import sys
sys.path.append("./scripts")
sys.path.append("./scripts/snakemake_helpers")

from snakemake_helpers import *

#----- Snakemake Set Up -----#

configfile: "hecatomb_config.yaml"

# Hecatomb DB paths
CONPATH = config["Paths"]["Contaminants"]
HOSTPATH = config["Paths"]["Host"]
BACPATH = config["Paths"]["Bacteria"]
AATARGET = config["Paths"]["TargetMMseqsAA"]
AATARGETCHECK = config["Paths"]["TargetMMseqsAACheck"]
NTTARGET = config["Paths"]["TargetMMseqsNT"]
NTTARGETCHECK = config["Paths"]["TargetMMseqsNTCheck"]

# Data paths
READDIR = config["Paths"]["Reads"]

# File paths
PHAGE = config["DatabaseFiles"]["Phage"]
NCBIACC = config["DatabaseFiles"]["NCBIAccession"] 

# Write out NCBI Accession path to be read in by R scripts
ncbiAccPath = open("./taxonomizr_ncbi_accession_path.txt", "w")
ncbiAccPath.write(NCBIACC)
ncbiAccPath.close()

# Java memory
XMX = config["System"]["Memory"]

# Tools
BBTOOLS = config["Tools"]["BBTools"]
R = config["Tools"]["R"]
SEQKIT = config["Tools"]["Seqkit"]
PULLSEQ = config["Tools"]["Pullseq"]
MMSEQS = config["Tools"]["MMseqs"]

#----- Rename input files if there are any to rename -----#

rename_files(config)

#----- Collect the Input Files -----#

# Pull sample names from the renamed R1 files in /data/renamed/ and store in a list
SAMPLES, = glob_wildcards(os.path.join(READDIR, "renamed", "{sample}_R1.fastq.gz"))

PATTERN_R1 = "{sample}_R1"
PATTERN_R2 = "{sample}_R2"

#----- Snakemake Workflow -----# 

#----- Contaminant Removal -----#
include: "contaminant_removal.snakefile"
#----- Cluster Count -----#
include: "cluster_count.snakefile"
#----- Merge Sequencing Tables -----#
include: "seqtable_merge.snakefile"
#----- MMSeqs2 Query Sequencing Table Against UniProt Viral Proteins Clustered at 99% (Virus Uniprot)-----#
# Extract phage (phage_tax_table.tsv), non-phage (viruses_seqs.fasta), and unclassified (pviral_aa_unclassified_seqs.fasta)
include: "mmseqs_pviral_aa.snakefile"
#----- MMSeqs2 Query Probable Viral Non-Phage Sequences (viruses_seqs.fasta) Against UniClust30 ProteinDB (Remove False Positives; Uni Plus Virus)-----#
include: "mmseqs_pviral_aa_check.snakefile"
#----- MMSeqs2 Query Probable Viral Unclassified Sequences (pviral_aa_unclassified_seqs.fasta) Against Refseq Virus NT UniVec Masked -----#
include: "mmseqs_pviral_nt.snakefile"
#----- Annotate Probable Viral Unclassified Sequences Search (Alignment) Results -----#
include: "mmseqs_pviral_nt_annotate.snakefile"
#----- MMSeqs2 Query Probable Non-Phage Viral Sequences Against UniClust30 + Virus UniProt ProteinDB (Remove False Positives) -----#
include: "mmseqs_pviral_nt_check.snakefile"
#----- Annotate Probable Viral Non-Phage Sequences Search (Alignment) Results -----#
include: "mmseqs_pviral_nt_check_annotate.snakefile"
#----- Concatenate Results -----#
include: "concatenate_results.snakefile"

rule all:
	input:
		os.path.join("results", "mmseqs_aa_out", "phage_tax_table.tsv"),
		os.path.join("results", "mmseqs_aa_out", "aln.m8"),
		os.path.join("results", "mmseqs_aa_checked_out", "aln.m8"),
		os.path.join("results", "mmseqs_aa_checked_out", "taxonomyResult.tsv"),
		os.path.join("results", "mmseqs_aa_checked_out", "viruses_checked_aa_tax_table.tsv"),
		os.path.join("results", "mmseqs_aa_checked_out", "unclassified_checked_aa_seqs.fasta"),
		os.path.join("results", "mmseqs_nt_checked_out", "phage_nt_seqs.fasta"),
		os.path.join("results", "viruses_tax_table.tsv"),
		os.path.join("results", "phage_tax_table.tsv"),
		os.path.join("results", "aa.aln.m8"),
		os.path.join("results", "nt.aln.m8")

rule clean:
	shell:
		"rm -rf ./QC/ ./clumped/"
