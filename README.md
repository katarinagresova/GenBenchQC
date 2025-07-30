# Genomic Benchmarks QC: Automated Quality Control for Genomic Machine Learning Datasets

# Table of Contents
- [Installation](#installation)
- [Running as command line tools](#running-as-command-line-tools)
- [Running from Python](#running-from-python)
- [Description](#description)
  - [Examples of running](#examples-of-running)
- [Development](#development)
- [License](#license)

## Installation

Install Genomic Benchmarks QC using pip:

```bash
pip install git+https://github.com/katarinagresova/GenBenchQC.git
```

One of the modes of running GenBenchQC (`evaluate_split`) requires additional dependency for clustering - [cd-hit](https://www.bioinformatics.org/cd-hit/cd-hit-user-guide). If you plan to run this mode, please install cd-hit as well either through conda:
```bash
conda install -c bioconda cd-hit
```
or follow the official [installation instructions](https://github.com/weizhongli/cdhit/wiki/2.-Installation).

## Running as command line tools

### Sequence evaluator

```bash
evaluate_sequences -h
usage: evaluate_sequences [-h] --input INPUT --format {fasta,csv,tsv} 
                          [--out_folder OUT_FOLDER] 
                          [--sequence_column SEQUENCE_COLUMN [SEQUENCE_COLUMN ...]]
                          [--label_column LABEL_COLUMN] [--label LABEL] 
                          [--report_types REPORT_TYPES [REPORT_TYPES ...]] 
                          [--end_position END_POSITION]
                          [--plot_type PLOT_TYPE]
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
  --report_types {json,html} [{json,html} ...]
                        Types of reports to generate. Default: [html]
  --end_position END_POSITION
                        End position of the sequences to plot in the per position plots.
  --plot_type {violin,boxen}
                        Type of the plot to generate for per sequence nucleotide content.
                        For bigger datasets, "boxen" is recommended. Default: boxen.
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
                        [--seq_report_types SEQ_REPORT_TYPES [SEQ_REPORT_TYPES]]
                        [--end_position END_POSITION]
                        [--plot_type PLOT_TYPE]
                        [--log_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--log_file LOG_FILE]

A tool for evaluating sequence datasets.

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
  --report_types {json,html,simple} [{json,html,simple} ...]
                        Types of reports to generate. Default: [html, simple].
  --seq_report_types {json,html} [{json,html} ...]
                        Types of reports to generate for individual groups of sequences. Default: [].
  --end_position END_POSITION
                        End position of the sequences to consider in per position statistics.
  --plot_type {violin,boxen}
                        Type of plot to use for visualizations. For bigger datasets, "boxen" in recommended.
                        Default: boxen.
  --log_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Logging level, default to INFO.
  --log_file LOG_FILE   Path to the log file.
```

### Train-test split evaluator

```bash
evaluate_split -h
usage: evaluate_split [-h] --train_input TRAIN_INPUT [TRAIN_INPUT ...]
                                  --test_input TEST_INPUT [TEST_INPUT ...] 
                                  --format {fasta,csv,tsv}
                                 [--sequence_column SEQUENCE_COLUMN [SEQUENCE_COLUMN ...]] 
                                 [--out_folder OUT_FOLDER] [--report_types REPORT_TYPES [REPORT_TYPES ...]]
                                 [--log_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--log_file LOG_FILE]

Check data leakage in dataset train-test split.

options:
  -h, --help            show this help message and exit
  --train_input TRAIN_INPUT [TRAIN_INPUT ...]
                        Path to the dataset file with training data. Can be multiple files that will be evaluated as one dataset part.
  --test_input TEST_INPUT [TEST_INPUT ...]
                        Path to the dataset file with testing data. Can be multiple files that will be evaluated as one dataset part.
  --format {fasta,csv,tsv}
                        Format of the input files.
  --sequence_column SEQUENCE_COLUMN [SEQUENCE_COLUMN ...]
                        Name of the columns with sequences to analyze for datasets in CSV/TSV format. Either one column or list of columns.
  --out_folder OUT_FOLDER
                        Path to the output folder.
  --report_types {json,html} [{json,html} ...]
                        Types of reports to generate. Default: [html]
  --log_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Logging level, default to INFO.
  --log_file LOG_FILE   Path to the log file.
  --identity_threshold IDENTITY_THRESHOLD
                        Identity threshold for clustering. Default: 0.95
  --alignment_coverage ALIGNMENT_COVERAGE
                        Alignment coverage for clustering. Default: 0.8
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
  END_POSITION,
  PLOT_TYPE)
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
  SEQ_REPORT_TYPES,
  END_POSITION,
  PLOT_TYPE)
```

### Train-test split evaluator

```python
from genbenchQC import evaluate_split

evaluate_split.run(
  TRAIN_INPUT_PATHS_LIST, 
  TEST_INPUT_PATHS_LIST, 
  FILE_FORMAT, 
  OUT_FOLDER,
  SEQUENCE_COLUMN_LIST, 
  REPORT_TYPES,
  IDENTITY_THRESHOLD,
  ALIGNMENT_COVERAGE)
```

## Description

**GenBenchQC** can be run in two ways: via command line or from Python. Both ways require the same input parameters and perform the same logic.

The tool proved three modes of running while checking for different things in the datasets. **evaluate_sequences** is used for generating report for a dataset (or its part) as one part and it counts some basic statistics with plots.
**evaluate_dataset** is used for evaluating datasets containing multiple classes and comparing, if different features of the sequences are similar across the different classes.
Finally, **evaluate_split** is used for checking data leakage between train-test splits.

You can choose to run the tool while having different dataset formats:
- **FASTA**: The input is a FASTA file / list of FASTA files. One file needs to contain sequences of one class if running *evaluate_sequences* mode.
- **CSV/TSV**: The input is a CSV/TSV file, and you provide the name of the column containing sequences. You can have either:
  - **multiple files**, each one containing sequences from one class (similar as with FASTA input)
  - **one file** containing sequences from multiple classes. In this case, when running *evaluate_sequences* mode, you need to provide the name of the column containing class labels so the tool can split the dataset into parts. The label classes can then be inferred, or you can specify their list by yourself. The dataset will then be split into pieces containing sequences with corresponding labels and analysis will be performed similarly as with multiple files.

When having CSV/TSV input, you can also decide to provide multiple sequence columns to analyze. In this case, the analysis in modes *evaluate_sequences* and *evaluate_dataset* will be performed for each column separately and lastly for sequences made by concatenating sequences throughout all the columns. 
*evaluate_split* mode will run only the concatenated sequences.

### Examples of running:

If you want to run the tool on example datasets, please clone the repository first to be able to access the example data:

```bash
git clone https://github.com/katarinagresova/GenBenchQC.git
cd GenBenchQC
```

Then you can run the following commands to evaluate example datasets.

#### Dataset evaluation

```bash
evaluate_dataset --input example_datasets/G4_positives.fasta example_datasets/G4_negatives.fasta --format fasta --out_folder example_outputs/G4_dataset
```

```bash
evaluate_dataset --input example_datasets/miRNA_mRNA_pairs_dataset.tsv --format tsv --sequence_column gene --label_column label --out_folder example_outputs/miRNA_mRNA_dataset --log_level DEBUG
```

```python
from genbenchQC import evaluate_dataset

evaluate_dataset.run(['example_datasets/miRNA_mRNA_pairs_dataset.tsv'], 'tsv', 'example_outputs/miRNA_mRNA_dataset', ['gene', 'noncodingRNA'], 'label', ['0', '1'])
```

**Outputs**

Running the above commands will create an output folder `example_outputs` containing the results of the evaluation. 
The output will contain (in the brackets, you can find the specific names when executing with the G4 example dataset):
- *dataset_report_\*.csv (dataset_report_label_G4_positives_vs_G4_negatives.csv)* - a CSV file with the results of the evaluation between two classes of the dataset. 
It contains information if the two dataset classes are similar enough and thus the test passed or failed for each feature.
- *dataset_report_\*.html (dataset_report_label_G4_positives_vs_G4_negatives.html)* - a report with the results of the evaluation between two classes of the dataset.
It contains basic information about the dataset and plots visualizing individual features of the dataset.
Each plot can be toned with red color, meaning the specific feature was too different between the two dataset classes and should be examined.
- *\*_plots (dataset_report_label_G4_positives_vs_G4_negatives_plots)* - folder containing all plots from the report in *.png* format
- *dataset_report_\*_duplicates.txt* - present only if the dataset contains duplicate sequences between classes. It contains a list of all the duplicate sequences.

If you trigger also reports for individual classes of the dataset (by providing format for the `seq_report_types` parameter or by running directly our `evaluate_sequences` tool), the output folder will contain descriptive report of the selected dataset part.

#### Train-test split evaluation

```bash
evaluate_split --train_input example_datasets/enhancers_train.csv --test_input example_datasets/enhancers_test.csv --format csv --out_folder example_outputs/enhancers_dataset --sequence_column sequence
```

**Outputs**
- *train_test_check_\*_report.[html, json] (train-test_check_enhancers_train_vs_enhancers_test_report.html)* - report with clusters containing sequences from both train and test dataset parts, that could indicate data leakage.

## Development

If you want to help with the development of Genomic Benchmarks QC, you are more than welcome to join in!

For a guidance, have a look at [CONTRIBUTING.md](CONTRIBUTING.md)

## License

Genomic Benchmarks QC is MIT-style licensed, as found in the LICENSE file.