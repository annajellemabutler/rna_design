# RNA Sequence Optimization 

# to-do list 04/25
* Check that code still runs on lax pattern
* Read up on how other papers have generated single-stable RNA molecules 
* For sid meeting note
  * Mfold webserver outputs a reduced (hopefully representative) list of suboptimal structures. I tried mfold on command line, doesn't run correctly, so would have to do manually to find sequences. 
  * All other options only have webservers (Quikfold, UNAfold), except for RNAstructure, which has similar output to RNAfold
  * Essentially, RNAsubopt uses MFE prediction, which is not necessarily to say that all the structures outputted are plausible
  
## PROGRESS UPDATES 

As of 04/10, the script works great. The next step is to figure out which pattern will (a) minimize the free energy of the optimal structure and (b) ensure the suboptimal structures are much less likely to exist. 

*** 

This repository contains code designed to identify primary RNA sequences which meet particular criteria regarding their length, GC content, and predicted secondary structure(s). Its intended application was in the design of synthetic RNA molecules on which a novel enzymatic probe for high-throughput RNA structure determination could be tested. The project formed part of a five-week laboratory rotation I completed in Prof. Sid Dey's lab at the University of California, Santa Barbara. 

## Overview
1. An RNA sequence is randomly generated according to user-defined rules for length and GC content.
2. The minimum free energy secondary structure is predicted and compared to the user-defined secondary structure (`desired_secondary_structure`).
3. If the structure is a match, the dG of the optimal structure is compared with the dG of the nextmost suboptimal secondary structure.
4. If the absolute difference in dG is greater than a user-defined value (or if there is no suboptimal structure to compare to), the sequence is returned. Else, the loop begins again. 

## Dependencies
* ViennaRNA

## Pipeline

## 1. Generating random RNA sequences 
The function, `random_RNA_generator`, generates random RNA sequences of customizable composition. It takes three inputs: `length` (desired length of generated sequence),`required_nucleotides` (a string containing the nucleotides whose count you want to control, e.g., `'GC'`), and `required_ratio` (the proportion of the total nucleotides which should be taken from the `required_nucleotides` set). 

## 2. Predicting secondary structure
The function, `predict_structure`, runs `RNAfold` on the generated RNA sequence. The output is `dot_bracket_structure` (the minimum free energy structure in dot-bracket notation) and `mfe` (the minimum free energy value in kcal/mole). 

## 3. Comparing dG with the next most-optimal structure 
If the predicted MFE structure matches the user-defined pattern, the function, `predict_subopt`, runs the RNA sequence through `RNAsubopt` to obtain the structure (`dot_bracket_structure`) and dG (`dG`) of the most optimal structure of the suboptimal ensemble. 

Arguments are . Be careful with -e; the number of structures returned grows exponentially with sequence length and energy range, and we only need the top two lines. 


