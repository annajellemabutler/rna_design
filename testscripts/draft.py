
#!/usr/bin/python

# RNA Sequence Search

# This script generates random RNA sequences of the specified length and GC content and filters for those whose predicted minimum free energy (MFE) structure matches the provided template (e.g., a single stem-loop with a stem of at least 20 nucleotides). Of the remaining sequences, only those whose second-most optimal structure deviates sufficently from the MFE structure are kept. The script currently terminates when a single sequence passing both these filters is found. 

import subprocess
from random import choice, shuffle
import re
import pandas

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
    output = 
##### MAIN LOOP #####

# Define desired secondary structure (single stem loop)
pattern = r'^\.*[(.]*\.*[).]*\.*$'
pattern_stringent = r'^\.*\({15,}\.{3,15}\){15,}\.*$'

# Initialize variable for looping purposes 
pattern_found = False

# Continue looping until a match is found
while not pattern_found:
    # Get random RNA sequence
    sequence = random_rna_generator(68, 'GC', 0.5)

    # Run RNAfold on the sequence to obtain MFE structure and dG
    output = open("RNA_fold_run.txt", "w")
    rnafold = subprocess.Popen(['RNAfold'], stdout=output, stdin=subprocess.PIPE)
    rnafold.communicate(input=sequence.encode())[0]

    # Extract relevant information (sequence, structure, MFE)
    with open('RNA_fold_run.txt', 'r') as file:
        first_line  = file.readline().strip() # This reads the RNA sequence 
        second_line = file.readline().strip() 

    dot_bracket, mfe = second_line.split(' ', 1) # This extracts the structure and MFE 

    # Compare the desired structure (pattern) to the observed structure
    match = re.search(pattern_stringent, dot_bracket)
    if match:
        pattern_found = True # MUST MOVE AFTER TESTING
        print("Appropriate MFE sequence found:", first_line,"\nStructure:",dot_bracket, "\ndG:", mfe)
        print("\nChecking suboptimal structures...")

        # Run RNAsubopt on identified sequence to get suboptimal structures
        with open("RNA_subopt_run.txt","w") as subopt_output:
           rnasub = subprocess.Popen(['RNAsubopt', '-e', '3','-s'], stdout=subopt_output, stdin = subprocess.PIPE)
           rnasub.communicate(input=first_line.encode())[0]

        # Extract relevant info from subopt output (structures, dGs) 
        with open("RNA_subopt_run.txt","r") as file:
            lines = file.readlines()

        sequence = lines[0].strip()
        dot_brackets = []
        energies = []
        for line in lines[1:]:
            dot_bracket, energy  = line.strip().split()
            dot_brackets.append(dot_bracket)
            energies.append(float(energy))

        # Turn into a dataframe for easier manipulation
        data = {
            'Dot-Bracket Structure': dot_brackets,
            'Minimum Free Energy': energies
        }

        df = pandas.DataFrame(data)
        print("\nComparison with Suboptimal Structure:")
        print(df.head(2))

        # Check if secondmost negative dG structure is at least 2 kcal/mole away from MFE structure 
        if abs(df.iloc[0]['Minimum Free Energy']-df.iloc[1]['Minimum Free Energy']) > 2:
            print("Success!")
            pattern_found = True # end the loop NEED TO MOVE THIS
        else:
            print("Too small difference.")          
    else:
        print("\nSearching...\n")

