# RNA Sequence Optimization 

This repository contains the code for a project aimed at identifying primary RNA sequences which meet particular criteria required for the testing of a novel enzymatic probe for high-throughput RNA structure determination. The project formed part of a five-week laboratory rotation I completed in Prof. Sid Dey's lab at the University of California, Santa Barbara. 

Specifically, this workflow achieves the following:
1. Generates random RNA sequences with user-defined length and GC content
2. Employs UNAFold OR UFold or RNAFold to predict the secondary structures and associated minimum free energy values (ΔG) for each sequence.
3. Returns only the sequences which either have only one predicted secondary structure (and one associated ΔG) or have two or more possible structures with differences in ΔG > 2 kcal/mol. 

## Usage
1. Clone the repository to your machine.
2. Install any dependencies required.
3. Execute the main program with your user-defined input to generate and analyze sequences.
4. Review the output files containing the sequences meeting the desired criteria.

## Dependencies
* ViennaRNA

# Pipeline

## 1. Generating random RNA sequences 
The script, `random_RNA_generator.py`, generates random RNA sequences with customizable composition. 

It takes four inputs:
* `length` (length of generated sequence),
* `required_nucleotides` (a string containing the nucleotides whose count you want to control),
* `required_ratio` (the proportion of the total nucleotides which are from the `required_nucleotides` set), and
* `num_sequences` (number of unique sequences to generate).

The output is a set containing the specified number of unique RNA sequences. 

Example use:
`sequences = random_rna_generator(10, 'GC', 0.5, 5)
print(sequences)`

This generates 5 unique RNA sequences of length 10 and with 50% GC content.


## 2. RNAFold
`RNAsubopt` calculates all suboptimal secondary structures within a given energy range above the MFE structure. Be careful, because the number of structures returned grows exponentially with sequence length and energy range. 
* Need to determine an upper limit for the energy range (i.e., if the optimal structure is -20, do we want to know about structures that are -10, or -5?)
* Need to somehow convert the output into readable information, so that we can extract structures where the difference in free energy is > 2. 
