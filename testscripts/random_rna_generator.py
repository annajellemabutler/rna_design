import subprocess
from random import choice, shuffle

def random_rna_generator(length, required_nucleotides, required_ratio):
    required_nucleotide_count = round(required_ratio * length) 
    complementary_nucleotides = list(set('GCUA') - set(required_nucleotides))
    sequence = [choice(required_nucleotides) for _ in range(required_nucleotide_count)]
    sequence += [choice(complementary_nucleotides) for _ in range(length - required_nucleotide_count)]
    shuffle(sequence)
    print(''.join(sequence))

