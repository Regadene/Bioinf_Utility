import os.path

from modules.dna_rna_tools import (
    transcribe,
    reverse,
    complement,
    reverse_complement,
)
from modules.fastq_utils import (
    is_seq_gc_in_bounds,
    is_seq_len_in_bounds,
    is_seq_quality_higher_than_threshold,
    read_fastq_seq_from_file,
    write_fastq_seq_to_file,
)


def run_dna_rna_tools(*seq_oper):
    """
    Applies a specified DNA/RNA operation to one or more sequences.

    This function takes a variable number of arguments, where the first but not
    last arguments are sequences (DNA/RNA), and the last argument is the
    operation to be applied to these sequences.
    Supported operations include 'transcribe', 'reverse', 'complement',
    and 'reverse_complement'.

    Parameters: *seq_oper: A variable-length argument list. - The first but
    not last arguments are expected to be sequences (str). - The last
    argument is the operation to apply (str). It must be one of the
    following: 'transcribe', 'reverse', 'complement', 'reverse_complement'.

    Returns: list or str: If a single sequence and operation are provided,
    the function returns the one string result. If multiple sequences are
    provided, an array of results for each sequence is returned. If the
    operation is invalid or the number of arguments is incorrect,
    the function return the error message.

    Example:
        result = run_dna_rna_tools("ATG", "transcribe")
        # Returns: "AUG"

        results = run_dna_rna_tools("ttG", "AT", "ATc", "complement")
        # Returns: ["aaC", "TA", "TAg"]
    """
    # definition of the supported operations dict
    operations = {
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }
    if len(seq_oper) == 2:
        sequence, operation = seq_oper[0], seq_oper[1]
        if operation in operations:
            return operations[operation](sequence)

    elif len(seq_oper) > 2:
        operation, sequences = seq_oper[-1], seq_oper[:-1]
        results = []
        for sequence in sequences:
            if operation in operations:
                operation_result = operations[operation](sequence)
                results.append(operation_result)
            else:
                return "Wrong operation"
        return results

    else:
        return "Wrong amount of arguments"


def filter_fastq(
    input_fastq,
    output_fastq=None,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0,
):
    """
    Filters FASTQ sequences based on GC content, sequence length, and quality
    score.

    This function reads sequences from an input FASTQ file and writes the
    filtered sequences to an output FASTQ file in a 'filtered' directory.
    Sequences are filtered according to specified GC content bounds, sequence
    length bounds, and a minimum quality threshold.

    Parameters:
    input_fastq (str):
        Path to the input FASTQ file containing the sequences to be filtered.

    output_fastq (str, optional):
        Name of the output file to which filtered sequences will be saved.
        If not specified, the output file will be saved with the same base name
        as the input file.
        File will be saved to the 'filtered' subfolder in FASTQ format.

    gc_bounds (tuple or int, optional):
        A tuple specifying the minimum and maximum GC-content percentages
        (default is (0, 100)). If an integer is provided, it is treated as the
        maximum GC content with a minimum of 0.

    length_bounds (tuple or int, optional):
        A tuple specifying the minimum and maximum sequence lengths
        (default is (0, 2**32)). If an integer is provided, it is treated as
        the maximum sequence length with a minimum of 0.

    quality_threshold (int, optional):
        The minimum acceptable quality score for sequences (default is 0).
        All sequences must have a quality score greater than or equal to this
        threshold to pass the filter.

    Returns:
        None: Writes the filtered sequences to the specified output FASTQ file.
        If any of the input parameters are invalid, returns the error message.

    Example:
        To filter sequences from a FASTQ file based on specific criteria:

        filter_fastq(
            "data/example_fastq.fastq",
            "each_filter_output.fastq",
            (40, 80),
            (10, 20),
            35
        )

        # This will write the filtered sequences to
        # 'filtered/each_filter_output.fastq'.
    """
    if type(gc_bounds) is int and gc_bounds > 0:
        gc_bounds = (0, gc_bounds)
    elif type(gc_bounds) is not tuple and len(gc_bounds) != 2:
        return "Wrong value for the gc_bounds argument"

    if type(length_bounds) is int and length_bounds > 0:
        length_bounds = (0, length_bounds)
    elif type(length_bounds) is not tuple and len(length_bounds) != 2:
        return "Wrong value for the length_bounds argument"

    if type(quality_threshold) is not int or quality_threshold < 0:
        return "Wrong value for the quality_threshold argument"

    if output_fastq is None:
        output_filtered_path = os.path.join(
            "filtered", os.path.basename(input_fastq)
        )
    else:
        output_filtered_path = os.path.join("filtered", output_fastq)

    if not os.path.exists("filtered"):
        os.mkdir("filtered")

    with (
        open(input_fastq, "r") as input_file,
        open(output_filtered_path, "w") as output_file,
    ):
        seq_name, sequence, quality_seq = read_fastq_seq_from_file(input_file)
        while seq_name:
            if (
                is_seq_gc_in_bounds(sequence, gc_bounds)
                and is_seq_len_in_bounds(sequence, length_bounds)
                and is_seq_quality_higher_than_threshold(
                    quality_seq, quality_threshold
                )
            ):
                write_fastq_seq_to_file(
                    seq_name, sequence, quality_seq, output_file
                )
            seq_name, sequence, quality_seq = read_fastq_seq_from_file(
                input_file
            )
