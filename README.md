# Bioinformatics Tools

## DNA/RNA Tools and FASTQ Filtering Script

This repository provides a Python script for working with DNA/RNA sequences and filtering FASTQ files based on various conditions like GC content, sequence length, and quality scores. It leverages external modules for sequence manipulation and filtering logic.

### Features

- **DNA/RNA Sequence Operations:**
  - Transcription from DNA to RNA.
  - Reversing sequences.
  - Computing complements of sequences.
  - Generating reverse complements.

- **FASTQ File Filtering:**
  - Reads from input FASTQ files and outputs filtered sequences to the file.
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

You can filter FASTQ sequences using the `filter_fastq()` function. This function can read sequences from an input FASTQ file and write the filtered sequences to an output FASTQ file.

```python
from bioinf_utility import filter_fastq

filter_fastq(
    "data/example_fastq.fastq",
    "each_filter_output.fastq",
    (40, 80),
    (10, 20),
    35
)

# This will write the filtered sequences to
# 'filtered/each_filter_output.fastq'.
```


### Notes

- Ensure your environment have the necessary permissions to create and write 
  subfolders and files during execution.
- Be noticed that output_fastq parameter of the filter_fastq is filename or 
  not specified, paths containing directories are not supported.