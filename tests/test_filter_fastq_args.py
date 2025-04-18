import pytest
from bioinf_utility import filter_fastq

DATA_DIR = "tests/_data"
FASTQ_FILE = DATA_DIR + "/example_fastq.fastq"

def test_invalid_gc_bounds_type():
    with pytest.raises(ValueError):
        filter_fastq(FASTQ_FILE, gc_bounds="bad")

def test_invalid_gc_bounds_list_order():
    with pytest.raises(ValueError):
        filter_fastq(FASTQ_FILE, gc_bounds=[60, 40])

def test_invalid_length_bounds_type():
    with pytest.raises(ValueError):
        filter_fastq(FASTQ_FILE, length_bounds="bad")

def test_invalid_length_bounds_list_order():
    with pytest.raises(ValueError):
        filter_fastq(FASTQ_FILE, length_bounds=[100, 50])

def test_invalid_quality_threshold():
    with pytest.raises(ValueError):
        filter_fastq(FASTQ_FILE, quality_threshold=-5)