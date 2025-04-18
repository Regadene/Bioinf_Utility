import os
import pytest
from bioinf_utility import filter_fastq

DATA_DIR = "tests/_data"
FASTQ_FILE = DATA_DIR + "/example_fastq.fastq"
OUTPUT_FILE_NAME = "test_out.fastq"
OUTPUT_FILE_PATH = f"filtered\{OUTPUT_FILE_NAME}"


def test_nonexistent_file():
    with pytest.raises(FileNotFoundError):
        filter_fastq(input_fastq="asdlkjhfsadlkjfyhwqpeuorfhypqwie.fastq")


def test_valid_file_filters_and_creates_output():
    filter_fastq(input_fastq=FASTQ_FILE, output_fastq=OUTPUT_FILE_NAME)

    assert os.path.exists(OUTPUT_FILE_PATH)
