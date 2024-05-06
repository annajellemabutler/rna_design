#!/usr/bin/python
# Working on main loop 

import subprocess
import sys
from random import choice, shuffle
import re

# Define sequence function
def random_rna_generator(length, required_nucleotides, required_ratio):
    required_nucleotide_count = round(required_ratio * length)
    complementary_nucleotides = list(set('GCUA') - set(required_nucleotides))
    sequence = [choice(required_nucleotides) for _ in range(required_nucleotide_count)]
    sequence += [choice(complementary_nucleotides) for _ in range(length - required_nucleotide_count)]
    shuffle(sequence)

    return ''.join(sequence)

# Define the pattern
pattern = r'^\.*[(.]*\.*[).]*\.*$'

# Initialize a variable to track whether the pattern is found
pattern_found = False

# Loop until the pattern is found
while not pattern_found:
    # Get random sequence
    sequence = random_rna_generator(68, 'GC', 0.5)

    # Run RNAfold on sequence
    output = open("RNA_fold_run.txt", "w")
    rnafold = subprocess.Popen(['RNAfold'], stdout=output, stdin=subprocess.PIPE)
    rnafold.communicate(input=sequence.encode())[0]
    output.close()

    # Open the file for reading
    with open('RNA_fold_run.txt', 'r') as file:
        MFE_sequence = file.readline().strip() 
        second_line = file.readline().strip()

    dot_bracket, mfe = second_line.split(' ', 1)

    # Check if the pattern is found
    match = re.search(pattern, dot_bracket)
    if match:
        pattern_found = True
        print("Pattern found. Sequence:", MFE_sequence,"Structure:",dot_bracket)
        # Run RNAsubopt on identified sequence
        subopt_output = open("RNA_subopt_run.txt","w")
        rnasub = subprocess.Popen(['RNAsubopt'], stdout=subopt_output, stdin = subprocess.PIPE)
        rnasub.communicate(input=MFE_sequence.encode())[0]
        output.close()
   
    else:
        print("Pattern not found. Structure:", dot_bracket)

