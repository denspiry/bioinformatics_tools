def is_valid_dna(seq):
    return all(base in 'ATCGatcg' for base in seq)


def is_valid_rna(seq):
    return all(base in 'AUCGaucg' for base in seq)


def is_nucleic_acid(seq):
    # Check if seq contains either DNA or RNA bases only
    has_t = any(base in 'Tt' for base in seq)
    has_u = any(base in 'Uu' for base in seq)
    return not (has_t and has_u) and (is_valid_dna(seq) or is_valid_rna(seq))


def reverse(seq):
    if not is_nucleic_acid(seq):
        return None
    return seq[::-1]


def complement(seq):
    if not is_nucleic_acid(seq):
        return None
    complement_map_dna = str.maketrans('ATCGatcg', 'TAGCtagc')
    complement_map_rna = str.maketrans('AUCGaucg', 'UAGCuagc')
    if is_valid_dna(seq):
        return seq.translate(complement_map_dna)
    else:
        return seq.translate(complement_map_rna)


def reverse_complement(seq):
    if not is_nucleic_acid(seq):
        return None
    return reverse(complement(seq))


def transcribe(seq):
    if not is_nucleic_acid(seq):
        return None
    transcription_map_dna = str.maketrans('ATCGatcg', 'AUCGaucg')
    if is_valid_dna(seq):
        # Replace T -> U
        return seq.translate(transcription_map_dna)
    else:
        return seq  # RNA is not transcribed


# Dictionary for matching procedures and functions
procedure_map = {
    'is_nucleic_acid': is_nucleic_acid,
    'transcribe': transcribe,
    'reverse': reverse,
    'complement': complement,
    'reverse_complement': reverse_complement
}
