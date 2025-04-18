import pytest
from biological_sequences import DNASequence, RNASequence, AminoAcidSequence

def test_dna_basic():
    dna = DNASequence("ATGC")
    assert len(dna) == 4
    assert dna[1].sequence == "T"
    assert dna[1:3].sequence == "TG"
    assert repr(dna) == "DNASequence(sequence=\"ATGC\")"
    assert dna.is_sequence_correct()

def test_dna_methods():
    dna = DNASequence("ATGC")
    assert dna.complement().sequence == "TACG"
    assert dna.reverse().sequence == "CGTA"
    assert dna.reverse_complement().sequence == "GCAT"

def test_dna_transcription():
    dna = DNASequence("ATTG")
    rna = dna.transcribe()
    assert isinstance(rna, RNASequence)
    assert rna.sequence == "AUUG"

def test_rna_sequence():
    rna = RNASequence("AUCG")
    assert len(rna) == 4
    assert rna[0:2].sequence == "AU"
    assert rna.is_sequence_correct()
    assert rna.complement().sequence == "TAGC"
    assert rna.reverse_complement().sequence == "CGAT"

def test_amino_acid_sequence():
    aa = AminoAcidSequence("MALWMRLLPL")
    assert aa.is_sequence_correct()
    assert len(aa) == 10
    assert aa.get_molecular_weight() == 1405.8

@pytest.mark.parametrize("seq, cls", [
    ("ATBGC", DNASequence),
    ("AUCB", RNASequence),
    ("XYZ", AminoAcidSequence),
])
def test_invalid_sequences(seq, cls):
    instance = cls(seq)
    assert not instance.is_sequence_correct()
