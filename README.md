# Detecting and Mitigating Biases in Genomic Data

## Installation

```bash
pip install git+https://github.com/katarinagresova/MitigatingBiases.git
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
from genData import evaluate_sequences

evaluate_sequences.run(INPUT_FASTA_PATH, OUT_FOLDER)
```

### Dataset evaluator

```python
from genData import evaluate_dataset

evaluate_dataset.run(POSITIVE_FASTA_PATH, NEGATIVE_FASTA_PATH, OUT_FOLDER)
```
