import os


def read_fastq_dict_from_file(input_file_path):
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
    at_header_line = input_file.readline().strip()
    sequence_line = input_file.readline().strip()
    input_file.readline().strip()
    quality_line = input_file.readline().strip()

    seq_name = at_header_line.strip().split(" ")[0]

    return seq_name, sequence_line, quality_line


def write_fastq_dict_to_file(sequences, output_fastq):
    with open(output_fastq, "w") as output_file:
        for seq_name, (sequence, quality_seq) in sequences.items():
            write_fastq_seq_to_file(
                seq_name, sequence, quality_seq, output_file
            )


def write_fastq_seq_to_file(seq_name, sequence, quality_seq, output_file):
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
    guanine_amount = seq.lower().count("g")
    cytosine_amount = seq.lower().count("c")
    gc_content = (guanine_amount + cytosine_amount) * 100 / len(seq)
    return gc_bounds[0] <= gc_content <= gc_bounds[1]


def is_seq_len_in_bounds(seq, length_bounds):
    return length_bounds[0] <= len(seq) <= length_bounds[1]


def is_seq_quality_higher_than_threshold(quality_seq, quality_threshold):
    total_quality_score = 0
    for character in quality_seq:
        total_quality_score += __get_phred_score(character)

    return total_quality_score / len(quality_seq) > quality_threshold


def __get_phred_score(character):
    return ord(character) - 33
