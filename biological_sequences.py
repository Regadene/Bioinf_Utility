from abc import ABC

# creating a dictionary for the complementary DNA sequence
COMPLEMENT_DICT_DNA = {
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
COMPLEMENT_DICT_RNA = {
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
COMPLEMENT_DICT_TRANS_DNA_TO_RNA = {
    "A": "A",
    "a": "a",
    "G": "G",
    "g": "g",
    "T": "U",
    "t": "u",
    "C": "C",
    "c": "c",
}

DNA_ALPHABET = set("ATGCatgc")
RNA_ALPHABET = set("AUGCaugc")
AMINO_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")

AMINO_MOLECULAR_WEIGHT_DICT = {
    "A":	89.1,
    "R":	174.2,
    "N":	132.1,
    "D":	133.1,
    "C":	121.2,
    "E":	147.1,
    "Q":	146.2,
    "G":	75.1,
    "H":	155.2,
    "I":	131.2,
    "L":	131.2,
    "K":	146.2,
    "M":	149.2,
    "F":	165.2,
    "P":	115.1,
    "S":	105.1,
    "T":	119.1,
    "W":	204.2,
    "Y":	181.2,
    "V":	117.1
}


class BiologicalSequence(ABC):
    alphabet: set[chr] = None

    def __init__(self, sequence: str):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.__class__(self.sequence[index])

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"{self.__class__.__name__}(sequence=\"{self.sequence}\")"

    def is_sequence_correct(self):
        return set(self.sequence).issubset(self.alphabet)


class NucleicAcidSequence(BiologicalSequence, ABC):
    complement_dict: dict = None

    def complement(self):
        if self.__class__ is NucleicAcidSequence:
            raise NotImplementedError("Use DNASequence or RNASequence instead.")
        return self.__class__("".join(self.complement_dict[n] for n in self.sequence))

    def reverse(self):
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    alphabet = DNA_ALPHABET
    complement_dict = COMPLEMENT_DICT_DNA

    def transcribe(self):
        return RNASequence("".join(COMPLEMENT_DICT_TRANS_DNA_TO_RNA[n] for n in self.sequence))


class RNASequence(NucleicAcidSequence):
    alphabet = RNA_ALPHABET
    complement_dict = COMPLEMENT_DICT_RNA


class AminoAcidSequence(BiologicalSequence):
    alphabet = AMINO_ALPHABET

    def get_molecular_weight(self):
        return round(sum(AMINO_MOLECULAR_WEIGHT_DICT[aa] for aa in self.sequence), 1)