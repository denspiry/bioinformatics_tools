import os


from typing import Generator, TextIO, Dict, Tuple, Union, List, Optional


from modules.module2 import (
    calculate_gc_content,
    calculate_average_quality,
    in_bounds
)


from modules.module1 import procedure_map


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0
) -> None:
    """
    Filters reads from an input FASTQ file.
    Writes passing reads to an output file.

    Filtering criteria:
    - GC content must be within gc_bounds
    - Sequence length must be within length_bounds
    - Average quality must be >= quality_threshold

    Arguments:
    input_fastq: str - path to the input FASTQ file
    output_fastq: str - name of the output file to write (will be saved in 'filtered/')
    seqs: dict - {name: (sequence, quality)}
    gc_bounds: tuple or float - GC content bounds in percent
    length_bounds: tuple or int - sequence length bounds
    quality_threshold: float - minimum average quality (Phred+33)
    """
    with open_output_file_safe(output_fastq) as (out_file, file_path):
        for name, sequence, plus, quality in read_fastq(input_fastq):
            gc_content = calculate_gc_content(sequence)
            avg_quality = calculate_average_quality(quality)
            seq_length = len(sequence)

        if (
            in_bounds(gc_content, gc_bounds)
            and in_bounds(seq_length, length_bounds)
            and avg_quality >= quality_threshold
        ):
            write_fastq_record(out_file, name, sequence, plus, quality)


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

    results: List[Optional[str]] = []

    for seq in sequences:
        if procedure_name == "is_nucleic_acid":
            results.append(is_nucleic_acid(seq))
            continue

        # Check if seq is nucleic acid
        if not is_nucleic_acid(seq):
            results.append(None)
            continue

        results.append(procedure_func(seq))

    return results[0] if len(results) == 1 else results