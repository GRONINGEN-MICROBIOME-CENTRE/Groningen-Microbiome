#!/bin/bash
#SBATCH --job-name=Diet_LLD
#SBATCH --error=DIET_err
#SBATCH --output=DIET_out
#SBATCH --mem=10gb
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1

python 01_Prepare_food_tables.py


