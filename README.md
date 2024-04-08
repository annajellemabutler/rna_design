# RNA Sequence Optimization 

## PROGRESS UPDATE 

3pm just worked out how to solo-out the dot-bracket structure (in tester_MFE.py). Now, want to incorporate that into find_MFE.py and see if that makes the pattern found. 
8am Currently trying to figure out the regex pattern to match dot-bracket structure indicating only one stem loop (working in `regex_pattern.py` on cluster). 

## TO-DO LIST
* download basic things to base environment (python, nano)
* 

This repository contains code designed to identify primary RNA sequences which meet particular criteria regarding their length, GC content, and predicted secondary structure(s). Its intended application was in the design of synthetic RNA molecules on which a novel enzymatic probe for high-throughput RNA structure determination could be tested. The project formed part of a five-week laboratory rotation I completed in Prof. Sid Dey's lab at the University of California, Santa Barbara. 

Specifically, this workflow achieves the following:
1. Generates random RNA sequences with user-defined length and GC content
2. Employs RNAFold to predict (a) the minimum free energy structure and (b) possible suboptimal secondary structures for each sequence. 
3. Filters for structures with only a single stem loop,  only the sequences which either have only one predicted secondary structure (and one associated ΔG) or have two or more possible structures with differences in ΔG > 2 kcal/mol. 

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
* Output of `RNAfold` is the MFE (most likely) secondary structure in dot-bracket notation.
* First step is to only continue the loop if the MFE contains only one stem loop. 
`RNAsubopt` calculates all suboptimal secondary structures within a given energy range above the MFE structure. Be careful, because the number of structures returned grows exponentially with sequence length and energy range. 
* Need to determine an upper limit for the energy range (i.e., if the optimal structure is -20, do we want to know about structures that are -10, or -5?)
* Need to somehow convert the output into readable information, so that we can extract structures where the difference in free energy is > 2.

***
* Finished Sun Apr 7 on python subopt.py (VSCODE) trying to figure out which regex pattern accurately describes what we want (single stem loop).
* Next steps I think is to turn the script into a loop, but using RNAfold instead of RNAsubopt. Stop the script when you find a sequence where the MFE structure only has one stem loop. Then get all the subopt structures and filter through those. 
