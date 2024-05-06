#!/usr/bin/python

# RNA Sequence Search

# Date: April 2024
# Project: Dey lab rotation

# Description
# This script generates random RNA sequences of the specified length and GC content and returns only those which match the following two conditions:
#       (1) The predicted minimum free energy (MFE) structure matches a user-defined pattern.
#       (2) The Gibbs free energy of the nextmost-optimal structure differs from that of the MFE structure by a user-defined value (or more).

# Main function input
#   - Secondary structure to search for in dot-bracket notation
#   - Required delta G (dG) difference between the optimal and nextmost-optimal predicted structures for a sequence
#   - Number of sequences to return

# Main function output
#   - Csv file containing matching sequences, their predicted structures and dGs, and the dGs of the nextmost-optimal structure
#   - The table is also printed to the terminal

# Imports
import subprocess
from random import choice, shuffle
import re
import pandas
import datetime

# Function which generates random RNA sequence of required length and GC content with poly-A tail
def random_rna_generator(length, required_nucleotides, required_ratio):
    adj_length = length - 6
    required_nucleotide_count = round(required_ratio * length)
    complementary_nucleotides = list(set("GCUA") - set(required_nucleotides))

    sequence = [choice(required_nucleotides) for _ in range(required_nucleotide_count)]
    sequence += [
        choice(complementary_nucleotides)
        for _ in range(length - required_nucleotide_count)
    ]
    shuffle(sequence)
    sequence += ["A"] * 6 
    return "".join(sequence)


# Function which runs RNAfold to predict MFE structure and extract relevant information for downstream use
def predict_structure(sequence):
    RNAfold_output = subprocess.check_output(
        ["RNAfold"], input=sequence.encode()
    ).decode()
    output_lines = RNAfold_output.strip().split("\n")
    dot_bracket_structure, mfe = output_lines[1].split(" ", 1)
    mfe = float(mfe.strip("()"))
    return dot_bracket_structure, mfe


# Function which runs RNAsubopt to get suboptimal secondary structures and extract relevant information for downstream use
def predict_subopt_structures(sequence):
    RNAsubopt_output = subprocess.check_output(
        ["RNAsubopt", "-e", "3", "-s"], input=sequence.encode()
    ).decode()
    output_lines = RNAsubopt_output.strip().split("\n")
    dot_bracket_structures = []
    gibbs_free_energies = []
    for line in output_lines[1:]:
        dot_bracket_structure, dG = line.strip().split()
        dot_bracket_structures.append(dot_bracket_structure)
        gibbs_free_energies.append(float(dG))
    return dot_bracket_structures, gibbs_free_energies


# Function which takes lists as input and compiles it into a dataframe which is printed and saved
def output_as_table(matching_sequences, matching_structures, matching_dGs, next_dGs):
    timestamp = datetime.datetime.now().strftime("%m-%d_%H-%M")  # Generate timestamp
    output_file = (
        f"/home/annajellema/gibbs/output/rna_sequence_set_{timestamp}.csv"  # Define output file name with timestamp
    )
    output = pandas.DataFrame(
        {
            "RNA Sequence": matching_sequences,
            "MFE Structure": matching_structures,
            "dG": matching_dGs,
            "Suboptimal dG": next_dGs,
        }
    )  # Can remove next_dGs; mainly included for testing purposes
    output.to_csv(output_file)
    print(f"\nOutput written to {output_file}\n")
    try:
        subprocess.run(["column", "-t", "-s", ",", output_file], check=True)
    except subprocess.CalledProcessError as e:  # Handle error
        print(f"Error: {e}")


# Main function; outputs RNA sequence(s) matching both criteria
def rna_sequence_search(
    desired_secondary_structure, percent_difference, required_sequences=1
):
    # Open empty lists
    matching_sequences = []
    matching_structures = []
    matching_dGs = []
    next_dGs = []

    # Main loop
    while len(matching_sequences) < required_sequences:
        rna_sequence = random_rna_generator(
            68, "GC", 0.5
        )  # Generate a random RNA sequence
        mfe_structure, mfe = predict_structure(
            rna_sequence
        )  # Predict its secondary structure
        match = re.search(
            desired_secondary_structure, mfe_structure
        )  # Check that it matches the desired secondary structure
        if match:
            suboptimal_structures, dGs = predict_subopt_structures(
                rna_sequence
            )  # Predict suboptimal secondary structures
            sorted_dGs = sorted(
                dGs
            )  # Sorting here rather than in RNAfold for efficiency
            # print("MFE structure:", suboptimal_structures[0], dGs[0])
            # print("Suboptimal structure:", suboptimal_structures[1], dGs[1])
            if (
                len(suboptimal_structures) == 1
            ):  # If there are no suboptimal structures, add sequence to list
                matching_sequences.append(rna_sequence)
                matching_structures.append(mfe_structure)
                matching_dGs.append(mfe)
                next_dGs.append("NA")
            elif (len(suboptimal_structures) > 1):
                threshold = abs(mfe) * percent_difference / 100
                if (abs(abs(mfe) - abs(dGs[1])) > threshold):  # If the dG difference between optimal and suboptimal structures is sufficient, add the sequence to the set
                    matching_sequences.append(rna_sequence)
                    matching_structures.append(mfe_structure)
                    matching_dGs.append(mfe)
                    next_dGs.append(sorted_dGs[1])
    output_as_table(matching_sequences, matching_structures, matching_dGs, next_dGs)

############################################################################################################

# Define desired secondary structure
# single stem loop, any structure: r"^\.*[(.]*\.*[).]*\.*$"  
# single stem loop, no bulges: r"^\.*\(*\.*\)*\.*$" 
# single stem of => 10 nt, any structure: r"^\.*[(.]{10,}\.*[).]{10,}\.*$" 
# single stem of => 8 nt, no bulges: r"^\.*\({8,}\.*\){8,}\.*$" 
# single or double stemloop: ^\.*\({3,}\.*\){3,}\.*(\({3,}\.*\){3,}\.*)?$



desired_secondary_structure = r"^\.*\({3,}\.*\){3,}\.*(\({3,}\.*\){3,}\.*)?$"

# Define percent difference
percent_difference = 20

# Define number of sequences
required_sequences = 5

# Call main function
rna_sequence_search(desired_secondary_structure, percent_difference, required_sequences)
