#!/bin/bash

#SBATCH --job-name=2_perfect_sequence
#SBATCH --time=48:00:00
#SBATCH --output=/home/annajellema/gibbs/output/2_perfect_rna_sequence_search_%j.out
#SBATCH --error=/home/annajellema/gibbs/output/2_perfect_rna_sequence_search_%j.err

cd /home/annajellema/gibbs
python3 rnadesign.py
