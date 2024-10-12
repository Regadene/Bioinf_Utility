# creating a dictionary for the complementary DNA sequence
complement_dict = {
    "A": "T",
    "a": "t",
    "G": "C",
    "g": "c",
    "T": "A",
    "t": "a",
    "C": "G",
    "c": "g",
}

# creating a dictionary for the complementary sequence
# from mRNA 5' - 3' to cDNA 3' - 5'
complement_dict_rna = {
    "A": "T",
    "a": "t",
    "G": "C",
    "g": "c",
    "U": "A",
    "u": "a",
    "C": "G",
    "c": "g",
}

# creating a dictionary for transcription DNA -> RNA (5' - 3' DNA ->
# 5' - 3' RNA)
complement_dict_trans_dna_to_rna = {
    "A": "A",
    "a": "a",
    "G": "G",
    "g": "g",
    "T": "U",
    "t": "u",
    "C": "C",
    "c": "c",
}


# check if the sequence is RNA
def is_rna(sequence):
    return set(sequence) <= {"A", "a", "U", "u", "G", "g", "C", "c"}


# check if the sequence is DNA
def is_dna(sequence):
    return set(sequence) <= {"A", "a", "T", "t", "G", "g", "C", "c"}


# transcription DNA -> RNA (transcription is performed only for DNA)
def transcribe(sequence):
    if is_dna(sequence):
        rna_sequence = ""
        for nucleotide in sequence:
            rna_sequence += complement_dict_trans_dna_to_rna[nucleotide]
        return rna_sequence


# reverses the DNA or RNA sequence
def reverse(sequence):
    return sequence[::-1]


# returns the complementary sequence
# (for DNA - complementary DNA, for RNA - cDNA)
def complement(sequence):
    if is_dna(sequence):
        complement_dna_sequence = ""
        for nucleotide in sequence:
            complement_dna_sequence += complement_dict[nucleotide]
        return complement_dna_sequence
    if is_rna(sequence):
        complement_rna_to_dna_sequence = ""
        for nucleotide in sequence:
            complement_rna_to_dna_sequence += complement_dict_rna[nucleotide]
        return complement_rna_to_dna_sequence


# returns the reverse complementary sequence
# (for DNA - DNA, for RNA - cDNA)
def reverse_complement(sequence):
    reversed_sequence = reverse(sequence)
    reversed_complement_sequence = complement(reversed_sequence)
    return reversed_complement_sequence
