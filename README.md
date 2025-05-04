# Genomic Benchmarks QC: Automated Quality Control for Genomic Machine Learning Datasets

## Installation

```bash
pip install git+https://github.com/katarinagresova/GenBenchQC.git
```

## Running as command line tools

### Sequence evaluator

```bash
evaluate_sequences INPUT_FASTA_PATH [--out_folder OUT_FOLDER]
```

### Dataset evaluator

```bash
evaluate_dataset POSITIVE_FASTA_PATH NEGATIVE_FASTA_PATH [--out_folder OUT_FOLDER]
```

## Running from Python

### Sequence evaluator

```python
from genbenchQC import evaluate_sequences

evaluate_sequences.run(INPUT_FASTA_PATH, OUT_FOLDER)
```

### Dataset evaluator

```python
from genbenchQC import evaluate_dataset

evaluate_dataset.run(POSITIVE_FASTA_PATH, NEGATIVE_FASTA_PATH, OUT_FOLDER)
```
