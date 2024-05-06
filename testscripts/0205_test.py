#!/usr/bin/python

##########################################################
# RNA Design

# Date: April 2024
# Project: Dey lab rotation, 

# Description:
# This script generates random RNA sequences of the specified length and GC content and returns only those which match the following two conditions:
#       (1) The predicted minimum free energy (MFE) structure matches a user-defined pattern.
#       (2) The Gibbs free energy of the nextmost-optimal structure differs from that of the MFE structure by at least the user-defined percentage.

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

def generate_rna_sequence():
    def generate_segment(length, gc_content):
        gc_count = int(length * gc_content)
        at_count = length - gc_count
        segment = ['G', 'C'] * gc_count + ['A', 'T'] * at_count
        random.shuffle(segment)
        return ''.join(segment[:length])

    # Define lengths and GC contents for each segment
    tail_5_len = 15
    stem_len = 15
    loop_len = 8
    tail_3_len = 2
    poly_a_len = 15
    stem_gc_content = 2/3
    overall_gc_content = 0.5

    # Generate initial segments
    stem_1 = generate_segment(stem_len, stem_gc_content)
    stem_2 = generate_segment(stem_len, stem_gc_content)
    loop = generate_segment(loop_len, overall_gc_content)
    tail_5 = generate_segment(tail_5_len, overall_gc_content)
    tail_3 = generate_segment(tail_3_len, overall_gc_content)
    poly_a = 'A' * poly_a_len

    # Ensure loop and tail_5 have at least 1 'C' and 2 'C's respectively
    if 'C' not in loop:
        loop = loop[:-1] + 'C'
    while tail_5.count('C') < 2:
        tail_5 = tail_5[:-1] + 'C'
    random.shuffle(list(tail_5))  # Shuffle to randomize the position of 'C'

    # Combine all parts
    rna_sequence = tail_5 + stem_1 + loop + stem_2 + tail_3 + poly_a

    # Adjust to meet overall GC content if necessary
    current_gc_count = rna_sequence.count('G') + rna_sequence.count('C')
    expected_gc_count = int(68 * overall_gc_content)
    while current_gc_count > expected_gc_count:
        # Replace GC with AT in non-stem areas if GC is too high
        rna_sequence = rna_sequence.replace('G', 'A', 1) if 'G' in rna_sequence else rna_sequence
        rna_sequence = rna_sequence.replace('C', 'T', 1) if 'C' in rna_sequence else rna_sequence
        current_gc_count = rna_sequence.count('G') + rna_sequence.count('C')
    while current_gc_count < expected_gc_count:
        # Replace AT with GC in non-stem areas if GC is too low
        rna_sequence = rna_sequence.replace('A', 'G', 1) if 'A' in rna_sequence else rna_sequence
        rna_sequence = rna_sequence.replace('T', 'C', 1) if 'T' in rna_sequence else rna_sequence
        current_gc_count = rna_sequence.count('G') + rna_sequence.count('C')

    return rna_sequence

# Generate an RNA sequence

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
def rna_sequence_search(desired_secondary_structure=r"^\.*\({3,}\.*\){3,}\.*$", dG_difference=1, required_sequences=1):
    matching_sequences = []
    matching_structures = []
    matching_dGs = []
    next_dGs = []
    count_all_sequences = 0
    count_sequences_with_matching_optimal = 0
    start_time = time.time()

    while len(matching_sequences) < required_sequences:
        rna_sequence = generate_rna_sequence() # generate sequence
        count_all_sequences += 1
        mfe_structure, mfe = predict_structure(
            rna_sequence
        )  # predict its secondary structure

        match = re.search(
            desired_secondary_structure, mfe_structure
        )  

        if match: # filter for those that match the desired secondary structure
            # Constraint: no G/C in stem base
            #base_index = mfe_structure.find('(')
            #base_nucleotides = rna_sequence[base_index:base_index+2]
            #if all(nucleotide not in ['G','C'] for nucleotide in base_nucleotides): 
            count_sequences_with_matching_optimal += 1
            suboptimal_structures, dGs = predict_subopt_structures(
                rna_sequence
            )  # predict suboptimal secondary structures
            sorted_dGs = sorted(
                dGs
            )  # sort here rather than in RNAfold for efficiency

            if (
                len(suboptimal_structures) == 1
            ):  # if there are no suboptimal structures (unlikely/impossible?) add sequence to list
                matching_sequences.append(rna_sequence)
                matching_structures.append(mfe_structure)
                matching_dGs.append(mfe)
                next_dGs.append("NA")
                count_sequences_passing_energy_requirements += 1

            elif (len(suboptimal_structures) > 1):
                print(f"Total count: {count_all_sequences}, Matching secondary: {count_sequences_with_matching_optimal}. \nLooking for structures above {sorted_dGs[0]}.")
                if abs(mfe - sorted_dGs[1]) > dG_difference:  # if the dG difference is sufficient, add the sequence to the set
                    matching_sequences.append(rna_sequence)
                    matching_structures.append(mfe_structure)
                    matching_dGs.append(mfe)
                    next_dGs.append(sorted_dGs[1])
                else:
                    print(f"dG too close: {sorted_dGs[1]}")

    output_as_table(matching_sequences, matching_structures, matching_dGs, next_dGs)

############################################################################################################

# Specify requirements
percent_difference = 0.2
required_sequences = 1
desired_secondary_structure = r"^\.*\({3,}\.*\){3,}\.*$"
    # single stem loop, any structure: r"^\.*[(.]*\.*[).]*\.*$"  
    # single stem loop, no bulges: r"^\.*\(*\.*\)*\.*$" 
    # single stem of => 10 nt, any structure: r"^\.*[(.]{10,}\.*[).]{10,}\.*$" 
    # single stem of => 8 nt, no bulges: r"^\.*\({8,}\.*\){8,}\.*$" 
    # single or double stemloop: r"^\.*\({3,}\.*\){3,}\.*(\({3,}\.*\){3,}\.*)?$"

# Call main function
rna_sequence_search(required_sequences=1, dG_difference=0.5)

# run a counter to see how many sequences you are running through

# notes from meeting
# if we know dG difference, 
# if 2 kcal/mole difference, then only 3% of molecules will in less stable structure (based on equilibrium constant)
# so lets order both -40 and -38, and -10 and -8 
# get all poly-Tail: 12-15 , so 53 plus 13 is 68 
# no. of Ts in cell-seq primer should guide number of As in poly-A

# do we want some with internal bulges/loops?
# make 5 or 6
# make 100 and then pull out 5 or 6 with different stem and loop lengths 
# specify how many Cs in the loop: make sure at least 1 C in the loop, but more the better (but at least 1)
# run it overnight to get lots, 1000
# look at bias in tail first rather than in loop, so no need to consider NC context
# at least 2 or 3 cytosines on each tail 
# remaining 5-6 in the stem 
# constraint: avoid having C in the first 2-3 bases of the stem (floppiness, = reaction is inefficient) 
# have 3 shorter stem, 3 longer stem (3 low stability and 3 high stability)
# but stem not shorter than 6-7 bases on each side 
# i.e., get different structures, with reasonable distribution of Cs
# order 6 molecules with that
#  POSSIBLy some other sequence on the 5' end (for TSO troubles): 10-12 nucleotides
#       - SOME Sprinkling of Ts, a few As and Gs, but not too many As
#       - leftover would be 41: so make stem 10 (each), loop of 7-8, and then 5-6 bases on right 
# cocnlusion was not to do the above, but might do in the future. # then 
# These are supposed to be DNA molecules ultimately 
# final crtieria
# - bias at AC (did not convert as well as we expected)
# - so make sure GC, AC, TC, etc. all combinations both in the single strand and in the stem to try and identify bias 
# - 