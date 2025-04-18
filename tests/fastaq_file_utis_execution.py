from bioinf_utility import filter_fastq

filter_fastq("data/example_fastq.fastq", "gc_bounds_output.fastq", (40, 80))
filter_fastq(
    "data/example_fastq.fastq",
    "length_bounds_output.fastq",
    length_bounds=(10, 20),
)
filter_fastq(
    "data/example_fastq.fastq", "quality_output.fastq", quality_threshold=35
)
filter_fastq(
    "data/example_fastq.fastq",
    "each_filter_output.fastq",
    (40, 80),
    (10, 20),
    35,
)
