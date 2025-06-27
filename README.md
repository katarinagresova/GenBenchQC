# Genomic Benchmarks QC: Automated Quality Control for Genomic Machine Learning Datasets

## Installation

Install Genomic Benchmarks QC using pip:

```bash
pip install git+https://github.com/katarinagresova/GenBenchQC.git
```

## Running as command line tools

### Sequence evaluator

```bash
evaluate_sequences -h
usage: evaluate_sequences [-h] --input INPUT --format {fasta,csv,tsv} 
                          [--out_folder OUT_FOLDER] 
                          [--sequence_column SEQUENCE_COLUMN [SEQUENCE_COLUMN ...]]
                          [--label_column LABEL_COLUMN] [--label LABEL] 
                          [--report_types REPORT_TYPES [REPORT_TYPES ...]] 
                          [--log_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}] 
                          [--log_file LOG_FILE]

A tools for evaluating sequence data.

options:
  -h, --help            show this help message and exit
  --input INPUT         Path to the input file.
  --format {fasta,csv,tsv}
                        Format of the input file.
  --sequence_column SEQUENCE_COLUMN [SEQUENCE_COLUMN ...]
                        Name of the columns with sequences to analyze for datasets in CSV/TSV format. 
                        Either one column or list of columns.
  --label_column LABEL_COLUMN
                        Name with the label column for datasets in CSV/TSV format.
                        Needed only if you want to select a specific class from the dataset.
  --label LABEL         Label of the class to select from the whole dataset. 
                        If not specified, the whole dataset is taken and analyzed as one piece.
  --out_folder OUT_FOLDER
                        Path to the output folder.
  --report_types REPORT_TYPES [REPORT_TYPES ...]
                        Types of reports to generate. Options: json, html. Default: [html]
  --log_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Logging level, default to INFO.
  --log_file LOG_FILE   Path to the log file.
```

### Dataset evaluator

```bash
evaluate_dataset -h
usage: evaluate_dataset [-h] --input INPUT [INPUT ...] --format {fasta,csv,tsv} 
                        [--out_folder OUT_FOLDER]
                        [--sequence_column SEQUENCE_COLUMN [SEQUENCE_COLUMN ...]] 
                        [--label_column LABEL_COLUMN] [--label_list LABEL_LIST [LABEL_LIST ...]]
                        [--regression]
                        [--report_types REPORT_TYPES [REPORT_TYPES ...]]
                        [--log_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--log_file LOG_FILE]

Evaluate positive and negative sequence datasets.

options:
  -h, --help            show this help message and exit
  --input INPUT [INPUT ...]
                        Path to the dataset file. 
                        Can be a list of files, each containing sequences from one class.
  --format {fasta,csv,tsv}
                        Format of the input files.
  --sequence_column SEQUENCE_COLUMN [SEQUENCE_COLUMN ...]
                        Name of the columns with sequences to analyze for datasets in CSV/TSV format. 
                        Either one column or list of columns.
  --label_column LABEL_COLUMN
                        Name with the label column for datasets in CSV/TSV format.
  --label_list LABEL_LIST [LABEL_LIST ...]
                        List of label classes to consider or "infer" to parse different 
                        labels automatically from label column. For datasets in CSV/TSV format.
  --regression          If True, label column is considered as a regression target and values
                        are split into 2 classes
  --out_folder OUT_FOLDER
                        Path to the output folder.
  --report_types REPORT_TYPES [REPORT_TYPES ...]
                        Types of reports to generate. Options: json, html, simple. 
                        Default: [html, simple].
  --log_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Logging level, default to INFO.
  --log_file LOG_FILE   Path to the log file.
```

## Running from Python

### Sequence evaluator

```python
from genbenchQC import evaluate_sequences

evaluate_sequences.run(
  INPUT_PATH, 
  FILE_FORMAT, 
  OUT_FOLDER, 
  SEQUENCE_COLUMN_LIST, 
  LABEL_COLUMN, 
  LABEL,
  REPORT_TYPES,
  LOG_LEVEL,
  LOG_FILE)
```

### Dataset evaluator

```python
from genbenchQC import evaluate_dataset

evaluate_dataset.run(
  INPUT_PATHS_LIST, 
  FILE_FORMAT, 
  OUT_FOLDER, 
  SEQUENCE_COLUMN_LIST, 
  LABEL_COLUMN, 
  LABEL_LIST,
  REGRESSION,
  REPORT_TYPES,
  LOG_LEVEL,
  LOG_FILE)
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

## Development

If you want to help with the development of Genomic Benchmarks QC, you are more than welcome to join in!

For a guidance, have a look at [CONTRIBUTING.md](CONTRIBUTING.md)

## License

Genomic Benchmarks QC is MIT-style licensed, as found in the LICENSE file.