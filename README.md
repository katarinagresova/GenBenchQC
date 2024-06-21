# Detecting and Mitigating Biases in Genomic Data

## Installation

```bash
pip install git+https://github.com/katarinagresova/MitigatingBiases.git@refactor
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

### Negative dataset generator

```bash
generate_negatives POSITIVE_FASTA_PATH NEGATIVE_FASTA_PATH
```
