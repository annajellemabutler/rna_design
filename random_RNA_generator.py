from random import choice, shuffle

def random_rna_generator(length, required_nucleotides, required_ratio, num_sequences):
    unique_sequences = set()
    while len(unique_sequences) < num_sequences:
        required_nucleotide_count = round(required_ratio * length) 
        complementary_nucleotides = list(set('GCUA') - set(required_nucleotides)) # create list containing non-required
        sequence = [choice(required_nucleotides) for _ in range(required_nucleotide_count)] # randomly select nucleotides from set 
        sequence += [choice(complementary_nucleotides) for _ in range(length - required_nucleotide_count)] # fill in the rest
        shuffle(sequence) # randomize the order
        unique_sequences.add(''.join(sequence))
    return unique_sequences
