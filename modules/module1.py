from typing import Dict, Tuple, Union, List, Optional


def is_valid_dna(seq: str) -> bool:
    """
    Checks if the sequence contains only DNA bases (A, T, C, G).

    Arguments:
    seq: str - nucleotide sequence

    Returns True if valid DNA, False otherwise.
    """

    return all(base in 'ATCGatcg' for base in seq)


def is_valid_rna(seq: str) -> bool:
    """
    Checks if the sequence contains only RNA bases (A, U, C, G).

    Arguments:
    seq: str - nucleotide sequence

    Returns True if valid RNA, False otherwise.
    """

    return all(base in 'AUCGaucg' for base in seq)


def is_nucleic_acid(seq: str) -> bool:
    # Check if seq contains either DNA or RNA bases only
    """
    Checks if the sequence is either DNA or RNA (but not both T and U).

    Arguments:
    seq: str - nucleotide sequence

    Returns True if valid DNA or RNA, False otherwise.
    """

    has_t: bool = any(base in 'Tt' for base in seq)
    has_u: bool = any(base in 'Uu' for base in seq)
    return not (has_t and has_u) and (is_valid_dna(seq) or is_valid_rna(seq))


def reverse(seq: str) -> Optional[str]: # Optional because of None
    """
    Reverses a DNA or RNA sequence.

    Arguments:
    seq: str - nucleotide sequence

    Returns reversed sequence as str.
    Returns None if the sequence is not valid DNA or RNA.
    """

    if not is_nucleic_acid(seq):
        return None
    return seq[::-1]


def complement(seq: str) -> Optional[str]: # Optional because of None
    """
    Computes the complement of a DNA or RNA sequence.

    Arguments:
    seq: str - nucleotide sequence

    Returns complemented sequence as str.
    Returns None if the sequence is not valid DNA or RNA.
    """

    if not is_nucleic_acid(seq):
        return None
    complement_map_dna: dict = str.maketrans('ATCGatcg', 'TAGCtagc')
    complement_map_rna: dict = str.maketrans('AUCGaucg', 'UAGCuagc')
    if is_valid_dna(seq):
        return seq.translate(complement_map_dna)
    else:
        return seq.translate(complement_map_rna)


def reverse_complement(seq: str) -> Optional[str]: # Optional because of None
    """
    Computes the reverse complement of a DNA or RNA sequence.

    Arguments:
    seq: str - nucleotide sequence

    Returns reverse-complemented sequence as str.
    Returns None if the sequence is not valid DNA or RNA.
    """

    if not is_nucleic_acid(seq):
        return None
    return reverse(complement(seq))


def transcribe(seq: str) -> Optional[str]: # Optional because of None
    """
    Transcribes DNA into RNA.

    Arguments:
    seq: str - nucleotide sequence

    Returns RNA sequence as str.
    Returns None if the sequense is not valid DNA or RNA
    """

    if not is_nucleic_acid(seq):
        return None
    transcription_map_dna: dict = str.maketrans('ATCGatcg', 'AUCGaucg')
    if is_valid_dna(seq):
        # Replace T -> U
        return seq.translate(transcription_map_dna)
    else:
        return seq  # RNA is not transcribed


# Dictionary for matching procedures and functions
procedure_map: Dict[str, callable] = {
    'is_nucleic_acid': is_nucleic_acid,
    'transcribe': transcribe,
    'reverse': reverse,
    'complement': complement,
    'reverse_complement': reverse_complement
}
