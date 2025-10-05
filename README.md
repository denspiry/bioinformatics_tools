### DNA/RNA Tools & FASTQ Filter

This Python package provides basic utilities for working with DNA/RNA sequences and filtering FASTQ-formatted reads based on GC content, sequence length, and quality scores.

## Structure
bioinformatics_tools/
├── main.py
├── modules/
│   ├── module1.py   # DNA/RNA tools
│   └── module2.py   # FASTQ filtering tools

## Features

* DNA/RNA utilities ('module1.py'):
    + 'reverse', 'complement', 'reverse_complement'
    + 'transcribe', 'is_nucleic_acid', etc.

* FASTQ filtering ('module2.py'):
    + GC content calculation
    + Quality score parsing (Phred+33)
    + Read filtering based on customizable bounds

* Dispatcher function:
    + 'run_dna_rna_tools(*args)' in 'main.py' allows you to apply operations like:
```python
run_dna_rna_tools('AUGC', 'reverse') # 'CGAU'
```

## How to Use

1. Run tools from 'main.py':

```bash 
python main.py
```

2. Use functions separately:

```python
from modules.module1 import reverse, transcribe
from modules.module2 import filter_fastq
```

3. Example filtering:

```python
filtered = filter_fastq(
    seqs=my_fastq_dict,
    gc_bounds=(40, 60),
    length_bounds=(50, 150),
    quality_threshold=30
)
```

4. Example operations with DNA and RNA

```python
run_dna_rna_tools('TTUU', 'is_nucleic_acid') # False
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```


## Requirements

* Python 3.6+
* No external libraries required (only standard library)

## Notes

* GC and quality thresholds are inclusive.
* Quality is interpreted using Phred+33 encoding.
* 'run_dna_rna_tools' accepts multiple sequences and returns transformed output(s).
* 'run_dna_rna_tools' accepts bot lower and upper case strings and doesn't change it


