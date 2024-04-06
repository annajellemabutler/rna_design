# RNA Sequence Optimization 

This repository contains the code for a project aimed at identifying primary RNA sequences which meet particular criteria required for the testing of a novel enzymatic probe for high-throughput RNA structure determination. The project formed part of a five-week laboratory rotation I completed in Prof. Sid Dey's lab at the University of California, Santa Barbara. 

Specifically, this workflow achieves the following:
1. Generates random RNA sequences with user-defined length and GC content
2. Employs UNAFold to predict the secondary structures and associated minimum free energy values (ΔG) for each sequence.
3. Returns only the sequences which either have only one predicted secondary structure (and one associated ΔG) or have two or more possible structures with differences in ΔG > 2 kcal/mol. 
