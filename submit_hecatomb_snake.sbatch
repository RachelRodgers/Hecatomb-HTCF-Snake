#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate /ref/mtblab/software/snakemake

snakemake --profile slurm_hecatomb --until seqtable_merge
