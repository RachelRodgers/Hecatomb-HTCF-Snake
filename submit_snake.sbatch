#!/bin/bash

ml snakemake/5.10.0-python-3.6.5

snakemake --profile slurm_hecatomb
