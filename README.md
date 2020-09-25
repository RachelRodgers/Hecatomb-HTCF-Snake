# Hecatomb Snakemake for HTCF

1. Login to HTCF.
2. Move to your /scratch/ directory.
3. Clone the repo:
```
git clone --recurse-submodules https://github.com/RachelRodgers/hecatomb_htcf_snake.git
```
4. Make a directory to hold the snakemake profile:
```
mkdir -p ~/.config/snakemake/slurm_hecatomb
```
5. Copy the cluster submit and profile files to the appropriate locations:
```
cd hecatomb_htcf_snake
cp slurm_config/config.yaml ~/.config/snakemake/slurm_hecatomb
cp slurm_submit/*.py ~/.config/snakemake
```
6. Create a directory to hold your raw sequencing reads and move your data to into that directory.
7. Edit the hecatomb_config.yaml file to point to your data directory (under Paths: Reads) and edit the Read1 Read2 and Extension patterns as needed (under Patterns:).
8. Submit in one of two ways:
	a. With sbatch script (preferred):
	```
	sbatch submit_hecatomb_snake.sbatch
	```
	b. Interactively (better for troubleshooting):
	```
	# start an interactive session:
	srun --mem=48G --cpus-per-task=8 -J hecatomb -p interactive --constraint=cpu_E52650 --pty /bin/bash -l
	
	# load snakemake:
	ml snakemake/5.10.0-python-3.6.5
	
	# dry run (prints steps and stops):
	snakemake --profile slurm_hecatomb -n
	
	# production run (run steps):
	snakemake --profile slurm_hecatomb
	```

