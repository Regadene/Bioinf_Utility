# Bioinformatics tools

## DNA/RNA Tools and FASTQ Filtering Script

This repository provides a Python script for working with DNA/RNA sequences and filtering FASTQ files based on various conditions like GC content, sequence length, and quality scores. It leverages external modules for sequence manipulation and filtering logic.

### Features

- **DNA/RNA Sequence Operations:**
  - Transcription from DNA to RNA.
  - Reversing sequences.
  - Computing complements of sequences.
  - Generating reverse complements.

- **FASTQ File Filtering:**
  - Filtering sequences based on GC content bounds.
  - Filtering by sequence length.
  - Filtering by quality score threshold.

### Usage

#### 1. DNA/RNA Sequence Operations

You can run various operations on DNA/RNA sequences using the `run_dna_rna_tools()` function. 

```python
from bioinf_utility import run_dna_rna_tools

# Example usage for a single sequence
result = run_dna_rna_tools("ATG", "transcribe")
# Returns: "AUG"

# Example usage for multiple sequences
results = run_dna_rna_tools("ttG", "AT", "ATc", "complement")
# Returns: ["aaC", "TA", "TAg"]
```

#### 2. FASTQ Filtering Script

You can filter fastq dictionaries using the `filter_fastq()` function.

```python
from bioinf_utility import filter_fastq

fastq_dict = {
    "@SRX079804": (
        "TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG",
        "FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB"
    ),
    "@SRX079810": (
        "CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC",
        "FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE"
    ),
}

# Applying the filter with specific bounds
filtered_dict = filter_fastq(
    fastq_dict, gc_bounds=(50, 65), length_bounds=(3, 78), quality_threshold=30
)

# Expected Output:
# {'@SRX079810': (
# 'CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC',
# 'FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE')}
```
