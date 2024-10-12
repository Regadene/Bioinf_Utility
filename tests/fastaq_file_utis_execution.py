from modules.fastq_utils import (
    write_fastq_dict_to_file,
    read_fastq_dict_from_file,
    convert_multiline_fasta_to_oneline,
    parse_blast_output,
)

from bioinf_utility import filter_fastq

from data.example_fastq_dict import EXAMPLE_FASTQ

fastq_dict = read_fastq_dict_from_file("data/example_fastq.fastq")
print(fastq_dict)
write_fastq_dict_to_file(fastq_dict, "output_data/output_fastq.fastq")
write_fastq_dict_to_file(EXAMPLE_FASTQ, "output_data/EXAMPLE_FASTQ.fastq")
# convert_multiline_fasta_to_oneline("data/example_multiline_fasta.fasta")
convert_multiline_fasta_to_oneline(
    "data/example_multiline_fasta.fasta",
    "output_data/onelinefasta_specified_path.fasta",
)
parse_blast_output(
    "data/example_blast_results.txt", "output_data/parsed_blast.fasta"
)


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
