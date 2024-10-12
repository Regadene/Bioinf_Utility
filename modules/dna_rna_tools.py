# создание словаря для комплементарной последовательности ДНК
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
# создание словаря для комплементарной последовательности
# из мРНК 5' - 3' в кДНК 3' - 5'
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
# создание словаря для транскрипции ДНК -> РНК (5' - 3' ДНК -> 5'-3' РНК)
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


# проверка на то, что последовательность - РНК
def is_rna(sequence):
    return set(sequence) <= {"A", "a", "U", "u", "G", "g", "C", "c"}


# проверка на то, что последовательность - ДНК
def is_dna(sequence):
    return set(sequence) <= {"A", "a", "T", "t", "G", "g", "C", "c"}


# транскрипция ДНК -> РНК (транскрипция осуществляется только для ДНК)
def transcribe(sequence):
    if is_dna(sequence):
        rna_sequence = ""
        for nucleotide in sequence:
            rna_sequence += complement_dict_trans_dna_to_rna[nucleotide]
        return rna_sequence


# разворачивает последовательность ДНК или РНК
def reverse(sequence):
    return sequence[::-1]


# возвращает комплементарную последовательность
# (для ДНК - комплементарную ДНК, для РНК - кДНК)
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


# возвращает обратную комплементарную последовательность
# (для ДНК - ДНК, для РНК - кДНК)
def reverse_complement(sequence):
    reversed_sequence = reverse(sequence)
    reversed_complement_sequence = complement(reversed_sequence)
    return reversed_complement_sequence
