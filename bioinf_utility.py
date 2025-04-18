import os.path
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import logging

from modules.dna_rna_tools import (
    transcribe,
    reverse,
    complement,
    reverse_complement,
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
    input_fastq: str,
    output_fastq: str = None,
    gc_bounds: int | list[int] = [0, 100],
    length_bounds: int | list[int] = [0, 2**32],
    quality_threshold: int = 0,
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

    gc_bounds (list or int, optional):
        A list specifying the minimum and maximum GC-content percentages
        (default is [0, 100]). If an integer is provided, it is treated as the
        maximum GC content with a minimum of 0.

    length_bounds (list or int, optional):
        A list specifying the minimum and maximum sequence lengths
        (default is [0, 2**32]). If an integer is provided, it is treated as
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
    logging.info("Checking input arguments")

    if not os.path.exists(input_fastq):
        raise FileNotFoundError("Input fastq file does not exists")

    if type(gc_bounds) is int and gc_bounds > 0:
        gc_bounds = (0, gc_bounds)
    elif type(gc_bounds) is not list \
            or len(gc_bounds) != 2 \
            or gc_bounds[0] > gc_bounds[1]:
        logging.error("Invalid gc_bounds: %s", gc_bounds)
        raise ValueError(
            "gc_bounds should be positive int or list of the 2 ints with the second value greater than first")

    if type(length_bounds) is int and length_bounds > 0:
        length_bounds = (0, length_bounds)
    elif type(length_bounds) is not list \
            or len(length_bounds) != 2 \
            or length_bounds[0] > length_bounds[1]:
        logging.error("Invalid length_bounds: %s", length_bounds)
        raise ValueError(
            "length_bounds should be positive int or list of the 2 ints with the second value greater than first")

    if type(quality_threshold) is not int or quality_threshold < 0:
        logging.error("Invalid quality_threshold: %s", quality_threshold)
        raise ValueError("quality_threshold should be positive int")

    if output_fastq is None:
        output_filtered_path = os.path.join(
            "filtered", os.path.basename(input_fastq)
        )
    else:
        output_filtered_path = os.path.join("filtered", output_fastq)

    if not os.path.exists("filtered"):
        os.mkdir("filtered")
        logging.info("Created 'filtered' directory.")
    logging.info("Input arguments are corrrect, start filtering process")

    filtered_records = []
    total_input_records_amount = 0
    for record in SeqIO.parse(input_fastq, "fastq"):
        total_input_records_amount += 1
        seq_len = len(record.seq)
        gc = gc_fraction(record.seq) * 100
        quality_seq = sum(record.letter_annotations["phred_quality"]) / seq_len

        if (
            length_bounds[0] <= seq_len <= length_bounds[1]
            and gc_bounds[0] <= gc <= gc_bounds[1]
            and quality_seq >= quality_threshold
        ):
            filtered_records.append(record)

    SeqIO.write(filtered_records, output_filtered_path, "fastq")
    logging.info(
        "Filtering complete. %d of %d sequences passed the filters and were saved to '%s'.",
        len(filtered_records), total_input_records_amount, output_filtered_path
    )
