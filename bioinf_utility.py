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
