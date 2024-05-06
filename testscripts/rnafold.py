#!/usr/bin/python
# trying to get MFE structure and filter based on that 

import subprocess
import sys
from random import choice, shuffle

# Define sequence function
def random_rna_generator(length, required_nucleotides, required_ratio):
    required_nucleotide_count = round(required_ratio * length)
    complementary_nucleotides = list(set('GCUA') - set(required_nucleotides))
    sequence = [choice(required_nucleotides) for _ in range(required_nucleotide_count)]
    sequence += [choice(complementary_nucleotides) for _ in range(length - required_nucleotide_count)]
    shuffle(sequence)

    return ''.join(sequence)

# Get random sequence 
sequence  = random_rna_generator(68,'GC', 0.5)

# Run RNAfold on sequence
output = open("RNA_fold_run.txt","w")
rnafold = subprocess.Popen(['RNAfold'],stdout=output, stdin=subprocess.PIPE)
rnafold.communicate(input=sequence.encode())[0]
output.close()

# pattern
pattern = r'\.*[(.]*\.*[).]*\.*'

# Open the file for reading
with open('RNA_fold_run.txt', 'r') as file:
    file_content = file.read()

if re.search(pattern, file_content):
    print("Pattern found.")
