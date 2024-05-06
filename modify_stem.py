import pandas as pd
import subprocess
import subprocess
from random import choice, shuffle
import random
import re
import pandas
import datetime
import time


def modify_sequence(sequence, dot_bracket):
    # Find loop and stem regions
    loop_end = dot_bracket.rfind('.', 0, dot_bracket.find(')')) + 1
    stem_start = dot_bracket.find('(') + 1
    print("Stem Start:", stem_start)
    stem_end = dot_bracket.rfind('(') + 1
    print("Stem End:", stem_end)
    # Define start of second part of stem following the loop
    second_stem_start = dot_bracket.find(')', loop_end) + 1
    print("Second Stem Start:", second_stem_start)
    second_stem_end = dot_bracket.rfind(')', loop_end) + 1
    print("Second Stem End:", second_stem_end)

    sequence_list = list(sequence)

    # Adjust first two nucleotides of the first stem to be A or U
    for nucleotide in range(stem_start - 1, stem_start + 1):
        if sequence_list[nucleotide] not in ['A', 'U']:
            sequence_list[nucleotide] = 'U' if sequence_list[nucleotide] == 'G' else 'A'

    # Adjust the last two nucleotides of the second stem to be A or U (complementary)
    for nucleotide in range(second_stem_end - 2, second_stem_end):
        if sequence_list[nucleotide] not in ['A', 'U']:
            sequence_list[nucleotide] = 'U' if sequence_list[nucleotide] == 'G' else 'A'

    return "".join(sequence_list)

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

# Function which compiles output lists into a dataframe which is printed and saved
def output_as_table(matching_sequences, matching_structures, matching_dGs, next_dGs):
    timestamp = datetime.datetime.now().strftime("%m-%d_%H-%M") 
    output_file = (
        f"/home/annajellema/gibbs/output/all_constraints_{timestamp}.csv"  
    )
    output = pandas.DataFrame(
        {
            "RNA Sequence": matching_sequences,
            "MFE Structure": matching_structures,
            "dG": matching_dGs,
            "Nextmost-optimal dG": next_dGs,
        }
    ) 
    output.to_csv(output_file)
    print(f"\nOutput written to {output_file}\n")
    subprocess.run(["column", "-t", "-s", ",", output_file], check=True)



# Main function which outputs RNA sequence(s) matching structure and energy criteria
def rna_final_constraints(input_file, dG_difference=2):
    
    df = pd.read_csv(input_file)

    matching_sequences = []
    matching_structures = []
    matching_dGs = []
    next_dGs = []

    for index, row in df.iterrows():
        modified_sequence = modify_sequence(row['RNA Sequence'], row['MFE Structure'])
        print(f"{row['RNA Sequence']}")
        print(f"{modified_sequence}")
        
    
        mfe_structure, mfe = predict_structure(modified_sequence)  # predict its secondary structure
        print(f"{mfe_structure}, {mfe}")

        suboptimal_structures, dGs = predict_subopt_structures(modified_sequence)  # predict suboptimal secondary structures
        sorted_dGs = sorted(dGs)

        if len(suboptimal_structures) == 1:  # if there are no suboptimal structures, add sequence to list
            matching_sequences.append(modified_sequence)
            matching_structures.append(mfe_structure)
            matching_dGs.append(mfe)
            next_dGs.append("NA")
        elif len(suboptimal_structures) > 1:
            if abs(mfe - sorted_dGs[1]) > dG_difference:  # if the dG difference is sufficient, add the sequence to the set
                matching_sequences.append(modified_sequence)
                matching_structures.append(mfe_structure)
                matching_dGs.append(mfe)
                next_dGs.append(sorted_dGs[1])
                print("SEQUENCE ADDED")
            else:
                print(f"dG too close: {sorted_dGs[1]}")

    output_as_table(matching_sequences, matching_structures, matching_dGs, next_dGs)

############################################################################################################

# Primary constraint check output file 
input_file = "/home/annajellema/gibbs/output/500_filtered_sequences_1.csv"
# Call main function
rna_final_constraints(input_file,dG_difference=2)

