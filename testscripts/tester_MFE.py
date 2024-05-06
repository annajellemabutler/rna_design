#!/usr/bin/python
# testing whether the reading second line part works

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
pattern = r'\.*[(.]*\.*[).]*\.*'

# Initialize a variable to track whether the pattern is found
pattern_found = False

# Get random sequence
sequence = random_rna_generator(68, 'GC', 0.5)

    # Run RNAfold on sequence
output = open("RNA_fold_run.txt", "w")
rnafold = subprocess.Popen(['RNAfold'], stdout=output, stdin=subprocess.PIPE)
rnafold.communicate(input=sequence.encode())[0]
output.close()
    
# Open the file for reading
with open('RNA_fold_run.txt', 'r') as file:
    next(file)
    second_line = file.readline().strip()

dot_bracket, mfe = second_line.split(' ',1)

print(dot_bracket)
