
#!/usr/bin/python

# RNA Sequence Search

# Testing rna_subopt output


import subprocess
from random import choice, shuffle
import re
import pandas
import datetime

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
    RNAsubopt_output = subprocess.check_output(['RNAsubopt','-e','3','-s'], input=sequence.encode()).decode() # Run RNAsubopt and store output as text
    output_lines = RNAsubopt_output.strip().split('\n')
    dot_bracket_structures = [] 
    gibbs_free_energies = []
    for line in output_lines[1:]:
        dot_bracket_structure, dG = line.strip().split()
        dot_bracket_structures.append(dot_bracket_structure)
        gibbs_free_energies.append(float(dG))
    return dot_bracket_structures, gibbs_free_energies

# Function to view output as a table for quick checking
def output_as_table(matching_sequences, matching_structures, matching_dGs, next_dGs):
    timestamp = datetime.datetime.now().strftime("%m-%d_%H-%M") # Generate timestamp
    output_file = f"/home/annajellema/gibbs/output/rna_sequence_set_{timestamp}.csv" # Define output file name with timestamp
    output = pandas.DataFrame({'RNA Sequence': matching_sequences, 'MFE Structure': matching_structures, 'dG': matching_dGs, 'Suboptimal dG': next_dGs})
    output.to_csv(output_file) 
    print(f"\nOutput written to {output_file}\n") # Print confirmation message
    try:
        subprocess.run(['column', '-t', '-s', ',', output_file], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")

# Main function to find RNA sequence matching all criteria 
def rna_sequence_search(desired_secondary_structure, dG_difference, required_sequences=1):
    matching_sequences = []
    matching_structures = []
    matching_dGs = []
    next_dGs = []
    while len(matching_sequences) < required_sequences:
        rna_sequence = random_rna_generator(68, 'GC', 0.5) # Generate a random RNA sequence 
        mfe_structure, mfe = predict_structure(rna_sequence) # Predict its secondary structure 
        match = re.search(desired_secondary_structure, mfe_structure) # Check that it matches the desired secondary structure
        if match:
            # print("\nMatching MFE structure found using sequence:", rna_sequence,".\nStructure:",mfe_structure, "\ndG:", mfe)
            # print("\nChecking suboptimal structures...")
            suboptimal_structures, dGs = predict_subopt_structures(rna_sequence) # Predict suboptimal secondary structures
            sorted_dGs = sorted(dGs) # sort here rather than in RNAfold for efficiency 
            # print("MFE structure:", suboptimal_structures[0], dGs[0])
            # print("Suboptimal structure:", suboptimal_structures[1], dGs[1])
            if len(suboptimal_structures) == 1: # If there are no suboptimal structures, add sequence
                matching_sequences.append(rna_sequence)
                matching_structures.append(mfe_structure)
                matching_dGs.append(mfe)
                next_dGs.append('NA')
                # print("Sucess! There are no suboptimal structures.")
            else: 
                if abs(sorted_dGs[0]-sorted_dGs[1]) > dG_difference: # If the dG difference between optimal and suboptimal structures is sufficient, add the sequence to the set 
                    matching_sequences.append(rna_sequence)
                    matching_structures.append(mfe_structure)
                    matching_dGs.append(mfe)
                    next_dGs.append(sorted_dGs[1])
                    #print("Success! The nextmost optimal structure deviates sufficiently from the MFE structure. Sequence added to set.")
    output_as_table(matching_sequences, matching_structures, matching_dGs, next_dGs)
    #return matching_sequences, matching_structures, matching_dGs

# Define desired secondary structure
stemloop_lax = r'^\.*[(.]*\.*[).]*\.*$' # matches single stem loop of any structure
stemloop_perfect = r'^\.*\(*\.*\)*\.*$' # matches single stem loop with no bulges or interior loops
stemloop_perfect_constrained = r'^\.*\({20,}\.*\){20,}\.*$' # matches single stem of ≥ 20 nucleotides with no bulges or loops 
stemloop_lax_constrained = r'^\.*[(.]{20,}\.*[).]{20,}\.*$' # matches single stem of ≥ 20 nucleotides including those with bulges and loops

# Call main function 
rna_sequence_search(stemloop_lax_constrained, 2, 3) # (pattern, required dG difference, number of sequences)
