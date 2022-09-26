#!/bin/bash

#eval "$(conda shell.bash hook)"
#conda activate snakemake

eval $( spack load --sh snakemake@3.11.2 )

snakemake --profile slurm_hecatomb --rerun-incomplete --until get_r1_singletons
