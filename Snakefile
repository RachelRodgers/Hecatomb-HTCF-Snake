"""
Master Snakefile for running the Hecatomb pipeline on HTCF.

Rachel Rodgers, Sep 2020
""" 

import os
import sys
sys.path.append("./scripts")

from hecatomb_helpers import *

#----- Snakemake Set Up -----#

configfile: "config.yaml"

# Hecatomb DB paths
CONPATH = config["Paths"]["Contaminants"]
HOSTPATH = config["Paths"]["Host"]
BACPATH = config["Paths"]["Bacteria"]
AATARGET = config["Paths"]["TargetMMseqsAA"]
NTTARGET = config["Paths"]["TargetMMseqsNT"]

# Data paths
READDIR = config["Paths"]["Reads"]

# File paths
PHAGE = config["DatabaseFiles"]["Phage"]

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
#----- MMseqs2 Query Viral Seqs Against AA DB -----#
include: "mmseqs_pviral_aa.snakefile"

rule all:
	input:
		os.path.join("results", "mmseqs_aa_out", "phage_tax_table.tsv"),
		os.path.join("results", "mmseqs_aa_out", "viruses_seqs.fasta"),
		os.path.join("results", "mmseqs_aa_out", "pviral_aa_unclassified_seqs.fasta")

rule clean:
	shell:
		"rm -rf ./QC/ ./clumped/"