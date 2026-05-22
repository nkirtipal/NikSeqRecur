# NikSeqRecur

NikSeqRecur is a recurrence-aware DNA sequence analysis framework designed for detecting approximate repeats, reverse-complement recurrences, and biologically meaningful mutation patterns in nucleotide sequences.

---

## Features

- Approximate repeat detection
- Reverse-complement recurrence analysis
- GC-content profiling
- Protein translation
- Synonymous / nonsynonymous mutation detection
- Stop gain / stop loss detection
- Transition / transversion classification
- Codon-position effect analysis
- Entropy-based filtering
- Genomic dotplot visualization
- Mutation overlap Venn diagrams
- BED export support

---

## Example Input Files

Example FASTA files are located in:

```text
example_data/
```

Included examples:

- `example_test.fa`
  - Synthetic DNA sequence for quick testing

- `sarscov2_spike.fasta`
  - SARS-CoV-2 spike gene example for biological recurrence analysis

---

## Example Outputs

Example outputs are located in:

```text
example_outputs/sarscov2/
```

Included outputs:

- Raw recurrence dotplots
- Filtered recurrence dotplots
- Mutation overlap Venn diagrams
- Recurrence distributions
- CSV exports
- BED genomic interval exports

---

## Installation

Install required packages:

```bash
pip install -r requirements.txt
```

Main dependencies:

- numpy
- pandas
- matplotlib
- biopython
- matplotlib-venn

---
## Tested Environment

NikSeqRecur was tested using:

- Python 3.13
- Linux environment

The current example workflow uses a 1000 bp test region for faster execution and visualization. Runtime and recurrence complexity may increase substantially for larger genomic regions depending on repeat density and analysis parameters.

Thresholds and filtering parameters can be adjusted according to the biological use case and computational requirements.

---
## Running NikSeqRecur

Run:

```bash
python nikseqrecur.py
```
---

## Output Types

NikSeqRecur generates:

- Recurrence result tables
- Mutation annotations
- Protein translations
- GC-content metrics
- Genomic visualizations
- BED interval files

---

## Project Goal

NikSeqRecur was developed as a flexible framework for recurrence-aware DNA sequence analysis and biologically informed mutation interpretation across diverse genomic datasets.

---

## Author

Nikhil K Kirtipal

GitHub:
https://github.com/nkirtipal
