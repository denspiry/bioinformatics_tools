import os
import re
from typing import Union, Set


def convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta: Union[str, None]
) -> None:
    """
    Converts a multi-line FASTA file into a single-line format for each sequence.

    Arguments:
    input_fasta: str - path to the input FASTA file
    output_fasta: str or None - path to the output file. 
    If not provided, creates <input_fasta>_oneline.fasta in the same directory.

    Returns:
    None. 

    Writes reformatted FASTA to the output file.
    """
    if output_fasta is None:
        base, ext = os.path.splitext(input_fasta)
        output_fasta = f"{base}_oneline.fasta"

    with open(input_fasta, 'r') as in_file, open(output_fasta, 'w') as out_file:
        sequence_lines = []
        header = None

        for line in in_file:
            line = line.strip()
            if not line:
                continue  # skip empty lines
            if line.startswith(">"):
                if header is not None:
                    # write previous header + sequence
                    out_file.write(header + '\n')
                    out_file.write(''.join(sequence_lines) + '\n')
                header = line # new header starts
                sequence_lines = [] # add sequence line in the new list
            else:
                sequence_lines.append(line)

        # Special part for the last record
        if header is not None:
            out_file.write(header + '\n')
            out_file.write(''.join(sequence_lines) + '\n')


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Parses a BLAST text output file and extracts the best hit (top match)
    for each query section.

    Arguments:
    input_file: str
        Path to the input BLAST result file (.txt).
    output_file: str
        Path to the output file where the extracted descriptions
        will be saved (one per line, sorted alphabetically).

    Returns:
    None

    Writes results to a text file. Each line contains the
    description of the best hit for one query.
    """
    # A set is used to remove duplicate descriptions
    results: Set[str] = set()

    with open(input_file) as file:
        lines = file.readlines()

    # This flag indicates that the next non-header line should be captured
    capture_next = False

    for line in lines:
        line = line.rstrip("\n")

        if "Sequences producing significant alignments:" in line:
            capture_next = True
            continue

        if capture_next:
            stripped = line.strip()

            # Skip blank lines and column headers
            if (
                not stripped
                or stripped.startswith("Scientific")
                or stripped.startswith("Description")
            ):
                continue

            # Split the line by two or more whitespace characters
            desc = re.split(r"\s{2,}", line.strip())[0]
            results.add(desc)
            capture_next = False
            continue

    sorted_results = sorted(results)

    with open(output_file, "w") as out_file:
        for name in sorted_results:
            out_file.write(name + "\n")