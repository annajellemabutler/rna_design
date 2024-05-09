#!/bin/bash

#SBATCH --job-name=_all6_sequences
#SBATCH --time=12:00:00
#SBATCH --output=/home/annajellema/gibbs/output/all6_slurm_rna_sequence_search_%j.out
#SBATCH --error=/home/annajellema/gibbs/output/all6_slurm_rna_sequence_search_%j.err

cd /home/annajellema/gibbs
python3 test.py
