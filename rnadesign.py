
#!/usr/bin/python

### Project: RNA Design (Dey lab rotation)
### Date: April 2024

# See README for description 

################################################################
### Imports & inputs ###

import subprocess
from random import choice, shuffle
import random
import re
import datetime
import pandas as pd 

# Define output directory 
output_directory = "/home/annajellema/gibbs/output"

# Define primary structure constraints
total_length = 68
poly_a_length = 13
min_stem_length = 12  # double-sided
min_cytosine_loop = 1
min_cytosine_tail = 2
min_cytosine_stem = 4//2 # single-sided)
min_loop_length = 6  

# Define energetic constraint (stemloop defined in loop)
dG_difference = 2

# Define number of sequences
batch_size=50 # number of sequences to generate before checking primary constraints
required_sequences = 1 # final number of sequences you want to get out

################################################################
### Functions ###

# Function which outlines the desired template, including biasing the stem towards stable stem loops
def rna_template(length, gc_bias=True):

    # This delineates regions of the molecule for increasing GC content in the "stem"
    polyA = poly_a_length
    eff_length = length - polyA
    min_tails = 7 
    max_tails = 10 
    tails_length = random.randint(min_tails, max_tails) # length of the unpaired nucleotides on either side of stem loop 
    stemloop_length = eff_length - (2 * tails_length) # length of stem loop 
    stem_length = (stemloop_length // 3)   
    loop_length = stemloop_length - (2 * stem_length)  # anything not tails_length and not stem is loop

    # Deal with GC content
    total_gc = length // 2 # number of GC needed to maintain overall 50% GC content
    if gc_bias: # make stem more GC-rich than rest of molecule for enhanced stability 
        stem_gc = (2 * stem_length * 2 // 3)  # make 2/3 of stem bases GC
    else:
        stem_gc = (2 * stem_length // 2)  # make 1/2 of the stem bases GC
    remaining_gc = total_gc - stem_gc # number of GC needed outside of stem 
    tail_gc = remaining_gc // 2  # evenly distribute remaining GC
    loop_gc = remaining_gc - tail_gc  

    # Generate random sequences for each subsection with specified GC content
    tail_5prime = generate_rna(tails_length+5, tail_gc // 2)
    tail_3prime = generate_rna(tails_length-5, tail_gc // 2)
    loop = generate_rna(loop_length, loop_gc)
    stem = generate_rna(stem_length, stem_gc // 2)
    complementary_stem = [complementary(nucleotide) for nucleotide in reversed(stem)]

    # Combine subsections 
    sequence = tail_5prime + stem + loop + complementary_stem + tail_3prime + (polyA *['A'])

    # Correct for length
    if len(sequence) < length:
        sequence = [random.choice(['A', 'C', 'G', 'U']) for _ in range(length-len(sequence)-1)] + sequence
    if len(sequence) > length:
        diff = len(sequence)-length
        sequence = sequence[diff:]

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

# Function which outputs a batch of sequences matching the secondary structure criteria
def secondary_struct_constraints(batch_size=50, dG_difference=2):
    desired_secondary_structure=r"^\.*\({3,}\.*\){3,}\.*$"

    matching_sequences = []
    matching_structures = []
    matching_dGs = []
    next_dGs = []

    count_all_sequences = 0
    count_stemloops = 0

    print(f"Looking for a batch of {batch_size} sequences which match secondary structure constraints...")

    while len(matching_sequences) < batch_size:
        rna_sequence = rna_template(total_length, gc_bias=True) # generate sequence
        count_all_sequences += 1

        mfe_structure, mfe = predict_structure(rna_sequence)  # predict its secondary structure
        match = re.search(desired_secondary_structure, mfe_structure) # check if stem loop

        if match: 
            count_stemloops += 1
            suboptimal_structures, dGs = predict_subopt_structures(rna_sequence)  # predict suboptimal secondary structures
            sorted_dGs = sorted(dGs)

            # check if energy constraint met 
            if abs(mfe - sorted_dGs[1]) > dG_difference:  # if the dG difference is sufficient, add the sequence to the set
                matching_sequences.append(rna_sequence)
                matching_structures.append(mfe_structure)
                matching_dGs.append(mfe)
                next_dGs.append(sorted_dGs[1])
            #else:
                #print(f"Energy difference insufficient: {sorted_dGs[1]}")
                
    #print(f"{batch_size} of the {count_stemloops} stemloop-forming sequences found (from {count_all_sequences} total sequences) match energy constraints. Now filtering batch for primary structure constraints...")
    
    return pd.DataFrame({
        'RNA Sequence': matching_sequences,
        'MFE Structure': matching_structures,
        'dG': matching_dGs,
        'Next dG': next_dGs
    })

# Function which checks batches of sequences for primary structure constraints
def primary_struct_constraints(results):
    valid_rows = []
    for index, row in results.iterrows():
        
        # Delineate regions
        dot_bracket = row["MFE Structure"]
        stem_start = dot_bracket.find('(') + 1
        stem_end = dot_bracket.rfind('(') + 1
        loop_start = dot_bracket.find('.', dot_bracket.rfind('(')) + 1
        loop_end = dot_bracket.rfind('.', 0, dot_bracket.find(')')) + 1
        stem_sequence = row["RNA Sequence"][stem_start-1:stem_end]
        loop_sequence = row["RNA Sequence"][loop_start-1:loop_end]
        tail_5 = row["RNA Sequence"][:stem_start-1]

        # Checking constraints (prints are mainly for debugging; can remove now)
        if len(row["RNA Sequence"]) > 68:
            print(f"Sequence {index} discarded: Incorrect length.")
            continue 

        if len(stem_sequence) < min_stem_length:
            print(f"Sequence {index} discarded: Stem length shorter than minimum required.")
            continue 
        
        if any(nucleotide in ['G', 'C'] for nucleotide in stem_sequence[:2]):
            print(f"Sequence {index} discarded: First two bases of the stem contains G or C.")
            continue

        #if stem_sequence[0] in ['G', 'C']:
            #print(f"Sequence {index} discarded: First base of the stem contains G or C.")
            #continue  # Uncomment to actually skip further processing of this sequence

        if stem_sequence.count('C') < min_cytosine_stem:
            print(f"Sequence {index} discarded: Insufficient cytosines in the stem.")
            continue

        if loop_sequence.count('C') < min_cytosine_loop: 
            print(f"Sequence {index} discarded: Insufficient cytosines in the loop.")
            continue

        if len(loop_sequence) < min_loop_length:
            print(f"Sequence {index} discarded: Insufficient nucleotides in the loop.")
            continue

        if tail_5.count('C') < min_cytosine_tail:
            print(f"Sequence {index} discarded: Insufficient cytosines in the 5' tail.")
            continue

        valid_rows.append(row)
        print(f"Sequence {index} added!")

    return pd.DataFrame(valid_rows)

# Main function (runs the above until required number of sequences found)
def main(required_sequences=1):
    collected_sequences = pd.DataFrame()
    while len(collected_sequences) < required_sequences:
        results = secondary_struct_constraints(batch_size=batch_size, dG_difference=dG_difference) 
        valid_results = primary_struct_constraints(results)
        collected_sequences = pd.concat([collected_sequences, valid_results])
        if len(collected_sequences) >= required_sequences:
            break

    return collected_sequences

################################################################
### Function call and output ###

test_set = main(required_sequences=required_sequences)

test_set = test_set.sort_values(by="dG")
timestamp = datetime.datetime.now().strftime("%m-%d_%H-%M") 
test_set.to_csv(f"{output_directory}/test_set_{timestamp}.csv", index=False)
print(f"{len(test_set)} sequences match all criteria; saved to '{output_directory}'")
