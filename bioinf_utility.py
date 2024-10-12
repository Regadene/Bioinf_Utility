from modules.dna_rna_tools import *


def run_dna_rna_tools(*seq_oper):
    # создание словаря с операциями
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
