#!/usr/bin/python

##########################################################
# RNA Design

# Date: April 2024
# Project: Dey lab rotation 

# Description:
# This script generates a set of random RNA sequences which match a set of rigorous constraints:
#   - Primary structure constraints
#       - 68 bp overall length
#       - 50% overall GC content
#       - 3' tail of 15 x adenosines 
#       - Reasonable distribution of cytosines 
#           - ≥ 1 cytosine in (what will be) the loop region
#           - 2-3 cytosines in (what will be) each tail 
#           - Remaining 5-6 in the stem 
#           - No cytosines in the first 2-3 bp of (what will be) the stem region
#       - Minimum of 6-7 bases per stem region 
#       - Both the tail and stem regions must have all combinations of NC (i.e., GC, AC, UC, CC). Don't worry about loop yet. 
#   - Energy & folding constraints 
#       - Predicted minimum free energy (MFE) structure is a single stem-loop (with no bulges or interior loops)
#       - Predicted suboptimal structures have absolute Gibbs free energy values ≤ 2 kcal/mole smaller than the MFE structure

# Main function input:
#       - Secondary structure to search for in dot-bracket notation (e.g., single stem loop)
#       - Required dG percent difference between optimal and nextmost-optimal predicted structures for a sequence
#       - Number of sequences to return

# Main function output:
#       - Csv file containing matching sequences, their predicted structures and dGs, and the dGs of the nextmost-optimal structures
#       - The table is also printed to the terminal



################################################################

import subprocess
from random import choice, shuffle
import random
import re
import pandas
import datetime
import time

# Function which biases RNA sequences towards forming stable stem loops
def rna_design(length, gc_bias=True):
    # SET THINGS 
    length = length 
    polyA = 13
    tail_length = int(3/13.5 * length)
    stem_length = int(3/13.5 * length)
    loop_length = length - (tail_length * 2) + stem_length

    # GENERATING RANDOM SEQUENCES
    total_gc = length // 2 # number of GC needed to maintain overall 50% GC content (34)
    if gc_bias: 
        stem_gc = (2 * stem_length * 2 // 3)  # make 2/3 of stem bases GC (20)
    else:
        stem_gc = (2 * stem_length // 2)  # make 1/2 of the stem bases GC
    remaining_gc = total_gc - stem_gc # number of GC needed outside of stem 
    tail_gc = remaining_gc // 2  
    loop_gc = remaining_gc - tail_gc  

    # Generate random sequences for each subsection with specified GC content
    tail_5prime = generate_rna(tail_length, tail_gc * 13 // 15)
    tail_3prime = generate_rna(tail_length-13, tail_gc * 2 // 15)
    tail_3prime = tail_3prime + (polyA * ['A'])
    loop = generate_rna(loop_length, loop_gc)
    stem = generate_rna(stem_length-2, stem_gc // 2)
    stem = ['A','U'] + stem
    complementary_stem = [complementary(nucleotide) for nucleotide in reversed(stem)]

    # Combine subsections 
    sequence = tail_5prime + stem + loop + complementary_stem + tail_3prime
    return "".join(sequence)

# Function which generates a random RNA sequence with specified GC content
def generate_rna(length, gc_count):
    gc_bases = ['G', 'C'] * (gc_count // 2)
    au_bases = ['A', 'U'] * ((length - gc_count) // 2)
    sequence = gc_bases + au_bases
    if len(sequence) < length:
        sequence.append(choice(['G', 'C'] if len(gc_bases) < gc_count else ['A', 'U']))
    shuffle(sequence)
    return sequence

# Function which returns the complementary RNA nucleotide
def complementary(nucleotide):
    return {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}.get(nucleotide, nucleotide)

# Function which runs RNAfold to predict MFE structure and extracts relevant information for later use
def predict_structure(sequence):
    RNAfold_output = subprocess.check_output(
        ["RNAfold"], input=sequence.encode()
    ).decode()
    output_lines = RNAfold_output.strip().split("\n")
    dot_bracket_structure, mfe = output_lines[1].split(" ", 1)
    mfe = float(mfe.strip("()")) # RNAfold outputs mfe in brackets
    return dot_bracket_structure, mfe

# Function which runs RNAsubopt to predict suboptimal secondary structures and extracts relevant information
def predict_subopt_structures(sequence):
    RNAsubopt_output = subprocess.check_output(
        ["RNAsubopt", "-e", "4", "-s"], input=sequence.encode()
    ).decode()
    output_lines = RNAsubopt_output.strip().split("\n")
    dot_bracket_structures = []
    gibbs_free_energies = []
    for line in output_lines[1:]:
        dot_bracket_structure, dG = line.strip().split()
        dot_bracket_structures.append(dot_bracket_structure)
        gibbs_free_energies.append(float(dG))
    return dot_bracket_structures, gibbs_free_energies

# Function which compiles output lists into a dataframe which is printed and saved
def output_as_table(matching_sequences, matching_structures, matching_dGs, next_dGs):
    timestamp = datetime.datetime.now().strftime("%m-%d_%H-%M") 
    output_file = (
        f"/home/annajellema/gibbs/output/rna_sequence_set_{timestamp}.csv"  
    )
    output = pandas.DataFrame(
        {
            "RNA Sequence": matching_sequences,
            "MFE Structure": matching_structures,
            "dG": matching_dGs,
            "Nextmost-optimal dG": next_dGs,
        }
    ) 
    output.to_csv(output_file)
    print(f"\nOutput written to {output_file}\n")
    subprocess.run(["column", "-t", "-s", ",", output_file], check=True)

# Main function which outputs RNA sequence(s) matching structure and energy criteria
def rna_sequence_search(desired_secondary_structure=r"^\.*\({3,}\.*\){3,}\.*$", dG_difference=2, required_sequences=1):
    
    matching_sequences = []
    matching_structures = []
    matching_dGs = []
    next_dGs = []
    count_all_sequences = 0
    count_sequences_with_matching_optimal = 0

    while len(matching_sequences) < required_sequences:
        rna_sequence = rna_design(68, gc_bias=True) # generate sequence
        count_all_sequences += 1
        mfe_structure, mfe = predict_structure(rna_sequence)  # predict its secondary structure
        match = re.search(desired_secondary_structure, mfe_structure)

        # Implementing various nucleotide constraints here (easier to match with structure) 
        if match: 
            # Constraint: no G/C in stem base
            base_index = mfe_structure.find('(')
            base_nucleotides = rna_sequence[base_index:base_index+2]
            if all(nucleotide not in ['G','C'] for nucleotide in base_nucleotides): 
                # Constraint: ≥ 1 C in loop 
                loop_start = mfe_structure.find('.', base_index)
                loop_end = mfe_structure.find(')',loop_start)
                loop_end = mfe_structure.rfind('.',loop_start,loop_end) + 1
                loop_sequence = rna_sequence[loop_start:loop_end]
                loop_Cs = loop_sequence.count('C')
                if loop_Cs >= 1:
                    # Constraint: ≥ 2 Cs in 5' tail
                    tail_end = mfe_structure.find('(') - 1
                    tail_nucleotides = rna_sequence[:tail_end]
                    if tail_nucleotides.count('C') >= 2:
                        count_sequences_with_matching_optimal += 1
                        suboptimal_structures, dGs = predict_subopt_structures(rna_sequence)  # predict suboptimal secondary structures
                        sorted_dGs = sorted(dGs)

                        if len(suboptimal_structures) == 1:  # if there are no suboptimal structures, add sequence to list
                            matching_sequences.append(rna_sequence)
                            matching_structures.append(mfe_structure)
                            matching_dGs.append(mfe)
                            next_dGs.append("NA")
                        elif len(suboptimal_structures) > 1:
                            if abs(mfe - dGs[1]) > dG_difference:  # if the dG difference is sufficient, add the sequence to the set
                                matching_sequences.append(rna_sequence)
                                matching_structures.append(mfe_structure)
                                matching_dGs.append(mfe)
                                next_dGs.append(dGs[1])
                            else:
                                print(f"dG too close: {dGs[1]}")
                else:
                    print(f"Not enough loop Cs: {loop_Cs}")
            else:
                print("Base nucleotides of stem include G or C. Skipping.")
    output_as_table(matching_sequences, matching_structures, matching_dGs, next_dGs)

############################################################################################################

# Call main function
rna_sequence_search(dG_difference=0.5)
