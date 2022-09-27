#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate /opt/apps/labs/mtblab/software/snakemake

#eval $( spack load --sh snakemake@3.11.2 )

snakemake --profile slurm_hecatomb --until seqtable_merge
