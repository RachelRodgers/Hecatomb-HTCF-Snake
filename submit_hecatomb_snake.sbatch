#!/bin/bash

eval "$(conda shell.bash hook)"

conda activate snakemake

snakemake --profile slurm_hecatomb --rerun-incomplete --until get_r1_singletons
