#!/usr/bin/python

# RNA Sequence Search

# This script generates random RNA sequences of the specified length and GC content and filters for those whose predicted minimum free energy (MFE) structure matches the provided template (e.g., a single stem-loop with a stem of at least 20 nucleotides). Of the remaining sequences, only those whose second-most optimal structure deviates sufficently from the MFE structure are kept. The script currently terminates when a single sequence passing both these filters is found. 

import subprocess
from random import choice, shuffle
import re

# Function which generates random RNA sequence of required length and GC content
def random_rna_generator(length, required_nucleotides, required_ratio):
    required_nucleotide_count = round(required_ratio * length) 
    complementary_nucleotides = list(set('GCUA') - set(required_nucleotides))
    sequence = [choice(required_nucleotides) for _ in range(required_nucleotide_count)]
    sequence += [choice(complementary_nucleotides) for _ in range(length - required_nucleotide_count)]
    shuffle(sequence)
    return ''.join(sequence)

# Function to run RNAfold to predict MFE secondary structure and extract relevant information for downstream use
def predict_structure(sequence):
    RNAfold_output = subprocess.check_output(['RNAfold'], input=sequence.encode()).decode() # Run RNAfold and store the output as text
    output_lines = RNAfold_output.strip().split('\n')
    dot_bracket_structure, mfe = output_lines[1].split(' ', 1)
    return dot_bracket_structure, mfe

# Function to run RNAsubopt to get suboptimal secondary structures and extract relevant information 
def predict_subopt_structures(sequence):
    RNAsubopt_output = subprocess.check_output(['RNAsubopt','e','1','-s'], input=sequence.encode()).decode() # Run RNAsubopt and store output as text
    output_lines = RNAsubopt_output.strip().split('\n')
    dot_bracket_structures = [] 
    gibbs_free_energies = []
    for line in output_lines[1:]:
        dot_bracket_structure, dG = line.strip().split()
        dot_bracket_structures.append(dot_bracket_structure)
        gibbs_free_energies.append(float(dG))
    return dot_bracket_structures, gibbs_free_energies

# Main function to find RNA sequence matching all criteria 
def rna_sequence_search(desired_secondary_structure, dG_difference, number_of_sequences_required=1):
    matching_sequences = set()
    while len(matching_sequences) < number_of_sequences_required:
        rna_sequence = random_rna_generator(68, 'GC', 0.5) # Generate a random RNA sequence 
        mfe_structure, mfe = predict_structure(rna_sequence) # Predict its secondary structure 
        match = re.search(desired_secondary_structure, mfe_structure) # Check that it matches the desired secondary structure
        if match:
            print("\nYay! Matching MFE structure found using sequence:", rna_sequence,".\nStructure:",mfe_structure, "\ndG:", mfe)
            print("\nChecking suboptimal structures...")
            suboptimal_structures, dGs = predict_subopt_structures(rna_sequence) # Predict suboptimal secondary structures
            print("MFE structure:", suboptimal_structures[0], dGs[0])
            print("Suboptimal structure:", suboptimal_structures[1], dGs[1])
            if abs(dGs[0]-dGs[1]) > dG_difference: # If the dG difference between optimal and suboptimal is sufficient, add the sequence to the set 
                matching_sequences.add(rna_sequence)
                print("Success! The nextmost optimal structure deviates sufficiently from the MFE structure. Sequence added to set.")
            else:
                print("The difference in dG is too small. Starting again...")
        else:
            print("\nSearching...\n")
    return matching_sequences

# Define desired secondary structure (single stem loop)
pattern = r'^\.*[(.]*\.*[).]*\.*$'
desired_secondary_structure = r'^\.*\({15,}\.{3,15}\){15,}\.*$'

# Call main function
matching_sequences = rna_sequence_search(desired_secondary_structure=pattern, dG_difference=0.5, number_of_sequences_required=3)
print(matching_sequences)
