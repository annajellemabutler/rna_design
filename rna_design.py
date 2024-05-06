
#!/usr/bin/python

### Project: RNA Design (Dey lab rotation)
### Date: April 2024

# See README for description and usage

################################################################
# IMPORTS AND FUNCTION DEFINITIONS

import subprocess
from random import choice, shuffle
import random
import re
import datetime
import pandas as pd 

# Function which biases RNA sequences towards forming stable stem loops
def rna_design(length, gc_bias=True):
    polyA = 13
    eff_length = length - polyA

    # Specify dimensions of stemloop 
    min_tails = 7 
    max_tails = 10 
    tails_length = random.randint(min_tails, max_tails) # length of the unpaired nucleotides on either side of stem loop 
    stemloop_length = eff_length - (2 * tails_length) # length of stem loop 
    stem_length = (stemloop_length // 3)   # customize desired length of stem here   
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
    if len(sequence) < length:
        sequence = [random.choice(['A', 'C', 'G', 'U']) for _ in range(length-len(sequence)-1)] + sequence

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

# Function which outputs RNA sequence(s) matching structure and energy criteria
def rna_sequence_search(desired_secondary_structure=r"^\.*\({3,}\.*\){3,}\.*$", dG_difference=1, required_sequences=1):
    matching_sequences = []
    matching_structures = []
    matching_dGs = []
    next_dGs = []
    count_all_sequences = 0
    count_sequences_with_matching_optimal = 0

    while len(matching_sequences) < required_sequences:
        rna_sequence = rna_design(68, gc_bias=True) # generate sequence
        count_all_sequences += 1
        mfe_structure, mfe = predict_structure(rna_sequence)  # predict its secondary structure
        match = re.search(desired_secondary_structure, mfe_structure)
        if match: 
            count_sequences_with_matching_optimal += 1
            suboptimal_structures, dGs = predict_subopt_structures(rna_sequence)  # predict suboptimal secondary structures
            sorted_dGs = sorted(dGs)

            if len(suboptimal_structures) == 1:  # if there are no suboptimal structures, add sequence to list
                matching_sequences.append(rna_sequence)
                matching_structures.append(mfe_structure)
                matching_dGs.append(mfe)
                next_dGs.append("NA")
            elif len(suboptimal_structures) > 1:
                if abs(mfe - sorted_dGs[1]) > dG_difference:  # if the dG difference is sufficient, add the sequence to the set
                    matching_sequences.append(rna_sequence)
                    matching_structures.append(mfe_structure)
                    matching_dGs.append(mfe)
                    next_dGs.append(sorted_dGs[1])
                #else:
                    #print(f"dG too close: {sorted_dGs[1]}")

    print(f"{required_sequences} sequences matching secondary structure constraints found. Checking primary structure constraints...\n")
    return pd.DataFrame({
        'RNA Sequence': matching_sequences,
        'MFE Structure': matching_structures,
        'dG': matching_dGs,
        'Next dG': next_dGs
    })

def validate_primary_sequence(results):
    valid_rows = []
    secondary_matches = len(results)
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

        # Checking constraints 
        if len(row["RNA Sequence"]) > 70:
            print(f"Sequence {index} discarded: Incorrect length.")
            continue 

        if len(stem_sequence) < min_stem_length:
            print(f"Sequence {index} discarded: Stem length shorter than minimum required.")
            continue 
        
        #if any(nucleotide in ['G', 'C'] for nucleotide in stem_sequence[:1]):
            #print(f"Sequence {index} discarded: First base of the stem contains G or C.")
            #continue

        if stem_sequence.count('C') < min_cytosine_stem:
            print(f"Sequence {index} discarded: Insufficient cytosines in the stem.")
            continue

        if loop_sequence.count('C') < min_cytosine_loop or len(loop_sequence) < min_loop_length:
            print(f"Sequence {index} discarded: Insufficient cytosines or nucleotides in the loop.")
            continue

        if tail_5.count('C') < min_cytosine_tail:
            print(f"Sequence {index} discarded: Insufficient cytosines in the 5' tail.")
            continue

        valid_rows.append(row)

    if valid_rows:
        output = pd.DataFrame(valid_rows)
        output = output.sort_values(by="dG")
        timestamp = datetime.datetime.now().strftime("%m-%d_%H-%M") 
        output_file = (f"/home/annajellema/gibbs/output/rna_sequence_{timestamp}.csv")
        output.to_csv(output_file, index=False)
        print(f"{len(output)} of {secondary_matches} sequences match primary structure criteria; saved to '{output_file}'")
        #subprocess.run(["column", "-t", "-s", ",", output_file], check=True)

    else:
        print("No valid sequences found.")

################################################################
# INPUTS 
 
# Define primary structure constraints
total_length = 68
poly_a_length = 13
min_stem_length = 12  # Total, 6 base pairs
min_cytosine_loop = 1
min_cytosine_tail = 2
min_cytosine_stem = 4//2 # because we halved stem)
min_loop_length = 6  # Minimum number of nucleotides in the loop

# Define secondary structure constraints
desired_secondary_structure = r"^\.*\({3,}\.*\){3,}\.*$"
dG_difference = 2
required_sequences = 10

################################################################
# Function calls 

results = rna_sequence_search(desired_secondary_structure=desired_secondary_structure, dG_difference=dG_difference,required_sequences=required_sequences) # get dataframe of sequences with secondary structure match
validate_primary_sequence(results) # save dataframe of sequences which match secondary and primary constraints
