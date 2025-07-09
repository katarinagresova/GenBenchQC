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

If you want to run the tool on example datasets, please clone the repository first:

```bash
git clone https://github.com/katarinagresova/GenBenchQC.git
cd GenBenchQC
```

Then you can run the following commands to evaluate example datasets.

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

#### Outputs

Running the above commands will create an output folder `example_outputs` containing the results of the evaluation.
The output will contain:
- *dataset_report_\*.csv* - a CSV file with the results of the evaluation between two classes of the dataset. 
It contains information if the two dataset classes are similar enough and thus the test passed or failed for each feature.
- *dataset_report_\*.html* - a report with the results of the evaluation between two classes of the dataset.
It contains basic information about the dataset and plots visualizing individual features of the dataset.
Each plot can be toned with red color, meaning the specific feature was too different between the two dataset classes and should be examined.
- *folder* - containing all plots from the report in *.png* format
- *dataset_report_\*_duplicates.txt* - present only if the dataset contains duplicate sequences between classes. It contains a list of all the duplicate sequences.

If you trigger also reports for individual classes of the dataset (by providing format for the `seq_report_types` parameter or by running directly our `evaluate_sequences` tool), the output folder will contain descriptive report of the selected dataset part.  
