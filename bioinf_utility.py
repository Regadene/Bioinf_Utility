from modules.dna_rna_tools import (
    transcribe,
    reverse,
    complement,
    reverse_complement,
)
from modules.fastq_filter_conditions import (
    is_seq_gc_in_bounds,
    is_seq_len_in_bounds,
    is_seq_quality_higher_than_threshold,
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
    seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0
):
    """
    Filters FASTQ sequences based on GC content, sequence length, and quality
    score.

    This function takes a dictionary of FASTQ sequences and filters them
    according to specified GC content bounds, sequence length bounds,
    and a minimum quality threshold.

    Parameters:
    seqs (dict):
        A dictionary where each key is a sequence name (str) and the value is
        a tuple consisting of (sequence, quality scores), where:
            - sequence (str): The DNA/RNA sequence.
            - quality scores (str): A string in phred+33 encoding representing
              the quality scores of the sequence.

    gc_bounds (tuple or int, optional):
        A tuple specifying the minimum and maximum GC-content percentages
        (default is (0, 100)). If an integer is provided, it is treated as the
        maximum GC content with a minimum of 0.

    length_bounds (tuple or int, optional):
        A tuple specifying the minimum and maximum sequence lengths (default is
        (0, 2**32)). If an integer is provided, it is treated as the maximum
        sequence length with a minimum of 0.

    quality_threshold (int, optional):
        The minimum acceptable quality score for sequences (default is 0). All
        sequences must have a quality score greater than or equal to this
        threshold to pass the filter.

    Returns:
        dict: A dictionary of filtered sequences that meet all
        conditions. The structure is the same as the input dictionary,
        with sequence names as keys and tuples of (sequence, quality scores) as
        values.

        str: An error message if any of the input parameters are invalid.

    Example:
        fastq_data = {
            "@SRX079804": (
                "TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG",
                "FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB"
            ),
            "@SRX079810": (
                "CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC",
                "FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE"
            ),
        }

        filtered_data = filter_fastq( fastq_data, gc_bounds=(50, 65),
        length_bounds=(3, 78), quality_threshold=30 )

        # Returns:
        # {'@SRX079810': (
        # 'CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC',
        # 'FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE')}

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

    filtered_fastq = {}

    for seq_name, (seq, quality_seq) in seqs.items():
        if (
            is_seq_gc_in_bounds(seq, gc_bounds)
            and is_seq_len_in_bounds(seq, length_bounds)
            and is_seq_quality_higher_than_threshold(
                quality_seq, quality_threshold
            )
        ):
            filtered_fastq[seq_name] = (seq, quality_seq)

    return filtered_fastq
