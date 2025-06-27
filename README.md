# Genomic Benchmarks QC: Automated Quality Control for Genomic Machine Learning Datasets

## Installation

```bash
pip install git+https://github.com/katarinagresova/GenBenchQC.git
```

## Running as command line tools

### Sequence evaluator

```bash
evaluate_sequences --input INPUT_PATH --format FASTA/CSV/TSV \
[--out_folder OUT_FOLDER, --sequence_column SEQUENCE_COLUMN, --label_column LABEL_COLUMN, --label LABEL, --log_level LOG_LEVEL, --log_file LOG_FILE]
```

### Dataset evaluator

```bash
evaluate_dataset --input INPUT_PATHS --format FASTA/CSV/TSV \
[--out_folder OUT_FOLDER, --sequence_column SEQUENCE_COLUMN, --label_column LABEL_COLUMN, --label_list LABEL1 LABEL2 ..., --regression, --log_level LOG_LEVEL, --log_file LOG_FILE]
```

## Running from Python

### Sequence evaluator

```python
from genbenchQC import evaluate_sequences

evaluate_sequences.run(INPUT_PATH, FILE_FORMAT, OUT_FOLDER, SEQUENCE_COLUMN_LIST, LABEL_COLUMN, LABEL)
```

### Dataset evaluator

```python
from genbenchQC import evaluate_dataset

evaluate_dataset.run(INPUT_PATHS_LIST, FILE_FORMAT, OUT_FOLDER, SEQUENCE_COLUMN_LIST, LABEL_COLUMN, LABEL_LIST, REGRESSION)
```

## Description

**GenBenchQC** can be run in two ways: via command line or from Python. Both ways require the same input parameters and perform the same logic.

You can choose to run the tool while having different dataset formats:
- **FASTA**: The input is a FASTA file / list of FASTA files. One file needs to contain sequences of one class.
- **CSV/TSV**: The input is a CSV/TSV file, and you provide the name of the column containing sequences. You can have either:
  - **multiple files**, each one containing sequences from one class (similar as with FASTA input)
  - **one file** containing sequences from multiple classes. In this case, you need to provide the name of the column containing class labels. The label classes can then be inferred, or you can specify their list by yourself. The dataset will then be split into pieces containing sequences with corresponding labels and analysis will be performed similarly as with multiple files.

When having CSV/TSV input, you can also decide to provide multiple sequence columns to analyze. In this case, the analysis will be performed for each column separately and lastly for sequences made by concatenating sequences throughout all the columns.

### Examples of running:

```bash
evaluate_dataset --input positives.fasta negatives.fasta --format fasta
```

```bash
evaluate_dataset --input dataset.csv --format csv --sequence_column seq --label_column label --labels 0 1 2
```

```python
from genbenchQC import evaluate_dataset

evaluate_dataset.run(['positives.tsv', 'negatives.tsv'], 'tsv', 'output_folder', ['seq1', 'seq2'])
```
