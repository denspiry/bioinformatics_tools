def calculate_gc_content(seq):
    """Calclulate GC content, return zero if empty string"""
    gc_count = sum(1 for base in seq if base in 'GCgc')
    return (gc_count / len(seq)) * 100 if seq else 0


def calculate_average_quality(quality_str):
    """Calculate read quality (phred33), return zero if empty string"""
    qualities = [ord(char) - 33 for char in quality_str]
    return sum(qualities) / len(qualities) if qualities else 0


def in_bounds(value, bounds):
    """Check if quality is in interval"""
    if isinstance(bounds, tuple):
        lower, upper = bounds
    else:
        lower, upper = 0, bounds
    return lower <= value <= upper
