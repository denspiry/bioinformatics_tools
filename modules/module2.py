from typing import Dict, Tuple, Union, List, Optional


def calculate_gc_content(seq: str) -> float:
    """
    Calculates the GC content of a nucleotide sequence as a percentage.

    Arguments:
    seq: str - nucleotide sequence

    Returns float - GC content in percent.
    Returns 0 if the sequence is empty.
    """

    gc_count = sum(1 for base in seq if base in 'GCgc')
    return (gc_count / len(seq)) * 100 if seq else 0


def calculate_average_quality(quality_str: str) -> float:
    """
    Calculates average read quality from a Phred+33 quality string.

    Arguments:
    quality_str: str - ASCII-encoded quality scores

    Returns float - average quality.
    Returns 0 if the string is empty.
    """

    qualities = [ord(char) - 33 for char in quality_str]
    return sum(qualities) / len(qualities) if qualities else 0


def in_bounds(value: Union[int, float], bounds: Union[Tuple[float, float], float]) -> bool:
    """
    Checks if a value is within given bounds (inclusive).

    Arguments:
    value: float or int
    bounds: tuple of (min, max) or single number as upper bound

    Returns True if value is within bounds, False otherwise.
    """

    if isinstance(bounds, tuple):
        lower, upper = bounds
    else:
        lower, upper = 0, bounds
    return lower <= value <= upper
