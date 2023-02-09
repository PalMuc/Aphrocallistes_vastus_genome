#!/bin/bash
#
#SBATCH --job-name=PINFISH
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --mem=40G
#SBATCH --qos=low
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL



snakemake --use-conda -j 40  all
