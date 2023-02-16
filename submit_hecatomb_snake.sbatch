#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate snakemake_7.22.0

snakemake --profile slurm_hecatomb
