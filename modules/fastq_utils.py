import os


def read_fastq_dict_from_file(input_file_path):
    """
    Read a FASTQ file and return a dictionary of sequences.
    """
    sequences_dict = {}

    with open(input_file_path, "r") as input_file:
        seq_name, sequence, quality_seq = read_fastq_seq_from_file(input_file)
        while seq_name:
            sequences_dict[seq_name] = (sequence, quality_seq)
            seq_name, sequence, quality_seq = read_fastq_seq_from_file(
                input_file
            )

    return sequences_dict


def read_fastq_seq_from_file(input_file):
    """
    Read a single sequence from a FASTQ file.
    """
    at_header_line = input_file.readline().strip()
    sequence_line = input_file.readline().strip()
    input_file.readline().strip()
    quality_line = input_file.readline().strip()

    seq_name = at_header_line.strip().split(" ")[0]

    return seq_name, sequence_line, quality_line


def write_fastq_dict_to_file(sequences, output_fastq):
    """
    Write a dictionary of sequences to a FASTQ file.
    """
    with open(output_fastq, "w") as output_file:
        for seq_name, (sequence, quality_seq) in sequences.items():
            write_fastq_seq_to_file(
                seq_name, sequence, quality_seq, output_file
            )


def write_fastq_seq_to_file(seq_name, sequence, quality_seq, output_file):
    """
    Write a single sequence to a FASTQ file.
    """
    output_file.writelines(
        [
            seq_name + "\n",
            sequence + "\n",
            "+" + seq_name[1:] + "\n",
            quality_seq + "\n",
        ]
    )


def convert_multiline_fasta_to_oneline(
    input_fasta_path, output_fasta_path=None
):
    """
    Convert a multi-line FASTA file to a single-line format.
    """
    if output_fasta_path is None:
        filename, extension = os.path.splitext(input_fasta_path)
        output_fasta_path = filename + "_oneline_format" + extension

    with (
        open(input_fasta_path, "r") as input_file,
        open(output_fasta_path, "w") as output_file,
    ):
        sequence_lines = []
        for line in input_file:
            line = line.strip()
            if line.startswith(">"):
                if sequence_lines:
                    output_file.write("".join(sequence_lines) + "\n")
                    sequence_lines = []
                output_file.write(line + "\n")
            else:
                sequence_lines.append(line)

        if sequence_lines:
            output_file.write("".join(sequence_lines) + "\n")


def parse_blast_output(input_file_path, output_file_path):
    """
    Parse a BLAST output file to extract best matches.
    """
    best_matches = set()

    with open(input_file_path, "r") as input_file:
        for line in input_file:
            line = line.strip()

            if line.startswith("Sequences producing significant alignments:"):
                input_file.readline()
                header_line = input_file.readline().strip()
                name_column_start_index = header_line.find("Name")
                best_match = (
                    input_file.readline()
                    .strip()[:name_column_start_index]
                    .strip()
                )
                best_matches.add(best_match)

    sorted_best_matches = sorted(best_matches)

    with open(output_file_path, "w") as output_file:
        for description in sorted_best_matches:
            output_file.write(description + "\n")


def is_seq_gc_in_bounds(seq, gc_bounds):
    """
    Check if the GC content of the sequence is within the specified bounds.
    """
    guanine_amount = seq.lower().count("g")
    cytosine_amount = seq.lower().count("c")
    gc_content = (guanine_amount + cytosine_amount) * 100 / len(seq)
    return gc_bounds[0] <= gc_content <= gc_bounds[1]


def is_seq_len_in_bounds(seq, length_bounds):
    """
    Check if the length of the sequence is within the specified bounds.
    """
    return length_bounds[0] <= len(seq) <= length_bounds[1]


def is_seq_quality_higher_than_threshold(quality_seq, quality_threshold):
    """
    Check if the average quality score of the sequence is above the threshold.
    """
    total_quality_score = 0
    for character in quality_seq:
        total_quality_score += __get_phred_score(character)

    return total_quality_score / len(quality_seq) > quality_threshold


def __get_phred_score(character):
    """
    Convert a quality score character to its Phred score.
    """
    return ord(character) - 33
