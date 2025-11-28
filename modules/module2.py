import os


from typing import Generator, TextIO, Dict, Tuple, Union, List, Optional


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


def read_fastq(file_path: str) -> Generator[Tuple[str, str, str, str], None, None]:
    """
    Reads a FASTQ file and yields one record (4 lines) at a time.

    Arguments:
    file_path: str - path to the input FASTQ file

    Yields:
    Tuple of (name, sequence, plus_line, quality) per read.
    """
    with open(file_path, 'r') as f:
        while True:
            lines = [f.readline().strip() for _ in range(4)]
            if not all(lines):
                break
            yield tuple(lines)


def write_fastq_record(
    output_file: TextIO,
    name: str,
    sequence: str,
    plus: str,
    quality: str
) -> None:
    """
    Writes a single FASTQ record (4 lines) to the output file.

    Arguments:
    output_file: TextIO - file to write to
    name: str - sequence identifier from fastq
    sequence: str - nucleotide sequence
    plus: str - '+' line
    quality: str - quality scores
    """
    output_file.write(f"{name}\n{sequence}\n{plus}\n{quality}\n")


def make_filtered_dir() -> str:
    """
    Creates a 'filtered' directory in the current working directory if it doesn't exist.

    Returns:
    str - path to the created or existing 'filtered' directory
    """
    filtered_dir = "filtered"
    os.makedirs(filtered_dir, exist_ok=True)
    return filtered_dir


def open_output_file_safe(output_filename: str) -> Tuple[TextIO, str]:
    """
    Opens an output file safely in the 'filtered' directory
    Doesn't overwrite existin file.

    Arguments:
    output_filename: str - name of the output file (not full path)

    Returns:
    Tuple[TextIO, str] - opened file and full path to the output file

    Raises:
    FileExistsError if the output file already exists.
    """
    filtered_dir = make_filtered_dir()
    full_path = os.path.join(filtered_dir, output_filename)

    if os.path.exists(full_path):
        raise FileExistsError(f"File '{full_path}' already exists.")

    return open(full_path, 'w'), full_path