import os
import shutil
import pytest
from Bio import SeqIO
from bioinf_utility import filter_fastq

DATA_DIR = "tests/_data"
FASTQ_FILE = DATA_DIR + "/example_fastq.fastq"
OUTPUT_FILE_NAME = "filtering_test_out.fastq"
OUTPUT_FILE_PATH = f"filtered\{OUTPUT_FILE_NAME}"


@pytest.mark.parametrize(
    "gc_bounds,length_bounds,quality_threshold,filtered_seqs_amount",
    [
        ([0, 100], [0, 10000], 0, 89),
        ([40, 60], [0, 1000], 0, 51),
        ([0, 100], [0, 50], 0, 50),
        ([0, 100], [50, 1000], 0, 40),
        ([0, 100], 10000, 35, 35),
        (60, 50, 25, 43),
    ],
)
def test_filter_fastq(
    gc_bounds, length_bounds, quality_threshold, filtered_seqs_amount
):
    filter_fastq(
        input_fastq=FASTQ_FILE,
        output_fastq=OUTPUT_FILE_NAME,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=quality_threshold,
    )

    assert os.path.exists(OUTPUT_FILE_PATH), "Output file was not created."

    records = list(SeqIO.parse(OUTPUT_FILE_PATH, "fastq"))
    assert (
        len(records) == filtered_seqs_amount
    ), f"Expected {filtered_seqs_amount} sequences, got {len(records)}"
