# Bioinformatics Tools

## DNA/RNA Tools and FASTQ Filtering Script

This repository provides a Python script for working with DNA/RNA sequences and filtering FASTQ files based on various conditions like GC content, sequence length, and quality scores. It leverages external modules for sequence manipulation and filtering logic.

### Features

- **DNA/RNA/Amino Sequence Classes:**
  - Transcription from DNA to RNA.
  - Reversing sequences.
  - Computing complements of sequences.
  - Generating reverse complements.
  - Counting molecular weight for Amino sequences

- **FASTQ File Filtering:**
  - Reads from input FASTQ files and outputs filtered sequences to the file.
  - Filtering by sequence length.
  - Filtering by quality score threshold.
  - CLI script support

#### 1. Biological Sequence Classes

This module provides object-oriented representations and utilities for handling DNA, RNA, and amino acid sequences. It includes methods for complementing, reversing, transcribing, and computing molecular weights of biological sequences.
##### Usage
```python
from biological_sequences import DNASequence, RNASequence, AminoAcidSequence

# Create a DNA sequence
dna = DNASequence("ATGCGT")

# Transcribe DNA to RNA
rna = dna.transcribe()

# Get reverse complement of DNA
rev_comp = dna.reverse_complement()
print(rev_comp) # Output: ACGCAT

# Create an RNA sequence and get its complement
rna_seq = RNASequence("AUGC")
print(rna_seq.complement()) # Output: TACG

# Create an amino acid sequence
protein = AminoAcidSequence("MALWMRLLPL")
weight = protein.get_molecular_weight()
print(weight)  # Output: 1405.8
```

#### 2. FASTQ Filtering Script

You can filter FASTQ sequences using the `filter_fastq()` function. This function can read sequences from an input FASTQ file and write the filtered sequences to an output FASTQ file.

```python
from bioinf_utility import filter_fastq

filter_fastq(
    "data/example_fastq.fastq",
    "each_filter_output.fastq",
    [40, 80],
    [10, 20],
    35
)

# This will write the filtered sequences to
# 'filtered/each_filter_output.fastq'.
```
#### 3. FASTQ Filtering CLI

The `filter_fastq_CLI.py` script allows you to filter FASTQ sequences directly from the command line. You can specify various options to control how the sequences are filtered based on GC content, sequence length, and quality score.

##### Basic Usage

To filter a FASTQ file, use the `-i` (or `--input`) argument to specify the input FASTQ file:

```bash
python filter_fastq_CLI.py -i tests/_data/example_fastq.fastq
```

You can request the CLI manual with the details about arguments:

```bash
python filter_fastq_CLI.py -h
```

##### Optional Parameters

- **GC Content Range**: Use the `--gc` argument to specify a GC content range (e.g., `--gc 30 60`).

```bash
python filter_fastq_CLI.py -i tests/_data/example_fastq.fastq --gc 30 60
```

- **Sequence Length Range**: Use the `-l` argument to filter sequences by their length (e.g., `-l 100 500`).

```bash
python filter_fastq_CLI.py -i tests/_data/example_fastq.fastq -l 100 500
```

- **Quality Score Threshold**: Use the `-q` argument to filter sequences by their quality score (e.g., `-q 30`).

```bash
python filter_fastq_CLI.py -i tests/_data/example_fastq.fastq -q 30
```

- **Output File**: Use the `--output` argument to specify the name of the output filtered FASTQ file.

```bash
python filter_fastq_CLI.py --input tests/_data/example_fastq.fastq --output filtered.fastq
```

##### Example Command

To filter sequences with GC content between 30 and 60, length between 100 and 500, and a quality score greater than or equal to 30, run:

```bash
python filter_fastq_CLI.py -i tests/_data/example_fastq.fastq --gc 20 60 -l 60 100 -q 30
```

##### Logging

All logs generated during the execution of the script, including any errors or warnings, are saved to a log file named `filter_fastq.log`. You can check this file for detailed information about the script's execution.


### Running Tests

To run the tests for the filtering logic and the CLI, follow these steps from root folder of the repository:

1. **Unit tests for the filtering function**: These tests ensure the `filter_fastq()` function works as expected by comparing the number of filtered sequences to the expected number.

    ```bash
    python -m pytest tests/test_filter_fastq_filtering.py tests/test_filter_fastq_files.py tests/test_filter_fastq_args.py
    ```

2. **CLI tests**: The `test_filter_fastq_CLI.py` module tests the command-line interface, ensuring that it correctly handles different arguments and returns the appropriate exit codes. It runs tests for both valid and invalid input arguments.

    ```bash
    python -m pytest tests/test_filter_fastq_CLI.py
    ```

3. **Run all tests**: To run both the unit tests and the CLI tests together:

    ```bash
    python -m pytest
    ```

This will execute both the filtering function and the CLI tests in a single command.


### Notes

- Ensure your environment have the necessary permissions to create and write 
  subfolders and files during execution.
- Be noticed that output_fastq parameter of the filter_fastq is filename or 
  not specified, paths containing directories are not supported.