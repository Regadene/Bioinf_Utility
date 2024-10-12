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
