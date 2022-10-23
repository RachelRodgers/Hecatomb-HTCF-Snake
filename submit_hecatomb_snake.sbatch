#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate snakemake_7.14.1

snakemake --profile slurm_hecatomb
