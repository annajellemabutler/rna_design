## testing chatgpt

import pandas as pd
import subprocess
import subprocess
from random import choice, shuffle
import random
import re
import pandas
import datetime
import time

# Load filtered sequences
df = pd.read_csv('/home/annajellema/gibbs/output/500_filtered_sequences_1.csv')

def modify_sequence(sequence, dot_bracket):
    # Find loop and stem regions
    loop_end = dot_bracket.rfind('.', 0, dot_bracket.find(')')) + 1
    stem_start = dot_bracket.find('(') + 1
    stem_end = dot_bracket.rfind('(') + 1
    # Define start of second part of stem following the loop
    second_stem_start = dot_bracket.find(')', loop_end) + 1
    second_stem_end = dot_bracket.rfind(')', loop_end) + 1

    sequence_list = list(sequence)
    print(f"Stem one: {stem_start}:{stem_end}, stem two: {second_stem_start}:{second_stem_end}")
    # Adjust first two nucleotides of the first stem to be A or U
    for nucleotide in range(stem_start - 1, stem_start + 1):
        if sequence_list[nucleotide] not in ['A', 'U']:
            sequence_list[nucleotide] = 'U' if sequence_list[nucleotide] == 'G' else 'A'

    # Adjust the last two nucleotides of the second stem to be A or U (complementary)
    for nucleotide in range(second_stem_end - 2, second_stem_end):
        if sequence_list[nucleotide] not in ['A', 'U']:
            sequence_list[nucleotide] = 'U' if sequence_list[nucleotide] == 'G' else 'A'

    print(sequence)
    print(dot_bracket)
    print("".join(sequence_list))

for index, row in df.iterrows():
        modified_sequence = modify_sequence(row['RNA Sequence'], row['MFE Structure'])


# Primary constraint check output file 
input_file = "/home/annajellema/gibbs/output/500_filtered_sequences_1.csv"
# Call main function
#rna_final_constraints(input_file,dG_difference=2)

