
#!/usr/bin/python

# RNA Sequence Search

# This script generates random RNA sequences of the specified length and GC content and filters for those whose predicted minimum free energy (MFE) structure matches the provided template (e.g., a single stem-loop with a stem of at least 20 nucleotides). Of the remaining sequences, only those whose second-most optimal structure deviates sufficently from the MFE structure are kept. The script currently terminates when a single sequence passing both these filters is found. 

import subprocess
from random import choice, shuffle
import re
import pandas
import time

# Function which generates random RNA sequence of required length and GC content
def random_rna_generator(length, required_nucleotides, required_ratio):
    required_nucleotide_count = round(required_ratio * length)
    complementary_nucleotides = list(set('GCUA') - set(required_nucleotides))
    sequence = [choice(required_nucleotides) for _ in range(required_nucleotide_count)]
    sequence += [choice(complementary_nucleotides) for _ in range(length - required_nucleotide_count)]
    shuffle(sequence)
    return ''.join(sequence)

# Function to run RNAfold and extract relevant information 
def run_rnafold(sequence):
    output = subprocess.check_output(['RNAfold'], input=sequence.encode()).decode()
    lines = output.strip().split('\n')
    sequence = lines[0].strip()
    dot_bracket, mfe = lines[1].split(' ', 1)
    return sequence, dot_bracket, mfe

# Function to run RNAsubopt and extract relevant information 
def run_subopt(sequence):
    output = subprocess.check_output(['RNAsubopt','e','1','-s'], input=sequence.encode()).decode()
    lines = output.strip().split('\n')
    sequence = lines[0].strip()
    dot_brackets = []
    energies = []
    for line in lines[1:]:
        dot_bracket, energy = line.strip().split()
        dot_brackets.append(dot_bracket)
        energies.append(float(energy))
    return sequence, dot_brackets, energies

# Main loop to find RNA sequence matching all criteria 
def rna_sequence_search(pattern_stringent, target_unique_sequences=1):
    unique_sequences = set()
    while len(unique_sequences) < target_unique_sequences:
        rna_sequence = random_rna_generator(68, 'GC', 0.5)
        sequence, dot_bracket, mfe = run_rnafold(rna_sequence)
        match = re.search(pattern_stringent, dot_bracket)
        if match:
            print("\nYay! Matching MFE structure found using sequence:", sequence,".\nStructure:",dot_bracket, "\ndG:", mfe)

            # Run RNA subopt 
            print("\nChecking suboptimal structures...")
            sequence, dot_brackets, energies = run_subopt(sequence)
            df = pandas.DataFrame({'Dot-Bracket Structure': dot_brackets, 'Minimum Free Energy': energies})
            print(df.head(2))

            # Compare free energies
            if abs(df.iloc[0]['Minimum Free Energy'] - df.iloc[1]['Minimum Free Energy']) > 0.1:
                print("Success! The nextmost optimal structure deviates sufficiently from the MFE structure.")
                print("\nSequence added to set.")
                unique_sequences.add(sequence)
            else:
                print("dG difference is too small. Starting again...")
        else:
            print("\nSearching...\n")
            time.sleep(1)

# Define desired secondary structure (single stem loop)
pattern = r'^\.*[(.]*\.*[).]*\.*$'
pattern_stringent = r'^\.*\({15,}\.{3,15}\){15,}\.*$'

# Call main function
rna_sequence_search(pattern)
   