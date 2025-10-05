from modules.module2 import (
    calculate_gc_content,
    calculate_average_quality,
    in_bounds
)


def filter_fastq(seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):
    """ Filter fastq-sequences"""
    filtered_seqs = {}

    for name, (sequence, quality) in seqs.items():
        gc_content = calculate_gc_content(sequence)
        avg_quality = calculate_average_quality(quality)
        seq_length = len(sequence)

        if (in_bounds(gc_content, gc_bounds) and
            in_bounds(seq_length, length_bounds) and
            avg_quality >= quality_threshold):
            filtered_seqs[name] = (sequence, quality)

    return filtered_seqs
