from typing import Dict, Tuple, Union, List, Optional


from modules.module2 import (
    calculate_gc_content,
    calculate_average_quality,
    in_bounds
)


from modules.module1 import procedure_map


def filter_fastq(
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0
) -> Dict[str, Tuple[str, str]]:
    """
    Filters FASTQ reads by GC content, length, and average quality.

    Arguments:
    seqs: dict - {name: (sequence, quality)}
    gc_bounds: tuple or float - GC content bounds in percent
    length_bounds: tuple or int - sequence length bounds
    quality_threshold: float - minimum average quality (Phred+33)

    Returns dict of filtered sequences with the same structure as input.
    """
    filtered_seqs = {}

    for name, (sequence, quality) in seqs.items():
        gc_content = calculate_gc_content(sequence)
        avg_quality = calculate_average_quality(quality)
        seq_length = len(sequence)

        if (
            in_bounds(gc_content, gc_bounds)
            and in_bounds(seq_length, length_bounds)
            and avg_quality >= quality_threshold
        ):
            filtered_seqs[name] = (sequence, quality)

    return filtered_seqs


def run_dna_rna_tools(*args: str) -> Union[str, List[Optional[str]]]:
    """
    Applies a selected nucleic acid operation to one or more sequences.

    Arguments:
    *args: one string or list of strings followed by a procedure name
    e.g., "ATCG", "complement"

    Returns:
    - Single result if one sequence is passed
    - List of results for multiple sequences
    Prints error messages for invalid input.
    """

    if len(args) < 2:
        print('ERROR: More arguments needed.')

    *sequences, procedure_name = args

    if procedure_name not in procedure_map:
        print(f"ERROR: Unknown procedure {procedure_name}")

    procedure_func = procedure_map[procedure_name]

    results: List[Optional[str]] = [procedure_func(seq) for seq in sequences]

    if len(results) == 1:
        return results[0] # to return one string
    return results
