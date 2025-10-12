import os


from typing import Union


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: Union[str, None]) -> None:
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