# Genomic Benchmarks QC: Automated Quality Control for Genomic Machine Learning Datasets

Genomic Benchmarks QC is a Python package and CLI toolkit for automated quality control of genomic datasets used in machine learning.
It helps detect biases, inconsistencies, and potential data leakage across sequences, dataset classes, and train-test splits — ensuring your datasets are reliable before model training.

## Features

### Provided Tools
- **evaluate_sequences** – QC of a single dataset or dataset subset.
- **evaluate_dataset** – QC across multiple dataset classes.
- **evaluate_split** – Train–test split leakage detection.

### General Features
- **Sequence-level QC** – Evaluate nucleotide composition, sequence length distribution, GC content, and more.
- **Class-level QC** – Compare multiple classes for feature similarity or bias.
- **Train–test split validation** – Detect potential data leakage through sequence similarity and clustering.
- **Multiple input format**s – Supports FASTA, CSV, and TSV datasets.
- **Customizable reporting** – Generate JSON, HTML, or simple text summaries.
- **Integration-ready** – Available as both CLI tools and a Python API.
- **Flexible sequence handling** – Works with single or multiple sequence columns.

## Installation

Install Genomic Benchmarks QC using pip:

```bash
pip install git+https://github.com/katarinagresova/GenBenchQC.git
```

If you plan to use `evaluate_split`, install [cd-hit](https://www.bioinformatics.org/cd-hit/cd-hit-user-guide):

```bash
conda install -c bioconda cd-hit
# or follow: https://github.com/weizhongli/cdhit/wiki/2.-Installation
```

## Quick Start

Clone the repository to access example datasets:

```bash
git clone https://github.com/katarinagresova/GenBenchQC.git
cd GenBenchQC
```

### Sequences QC

```bash
evaluate_sequences \
  --input example_datasets/G4_positives.fasta \
  --format fasta \
  --out_folder example_outputs/G4_dataset_positives
```

The same evaluation would be executed from Python as follows:

```python
from genbenchQC import evaluate_sequences

evaluate_sequences.run(
  input_file='example_datasets/G4_positives.fasta', 
  input_format='fasta',
  out_folder='example_outputs/G4_dataset_positives'
)
```

**Outputs**

Running the above commands for `evaluate_sequences` tool will create `example_outputs/G4_dataset_positives` folder with the following results:
- `G4_positives_plots` folder containing plots about different sequence features (nucleotide content, dinucleotide content, GC content, length, per-position nucleotide content).
- `G4_positives_report.html` file with descriptive statistics about sequences and including plots from `G4_positives_plots` folder.

### Dataset QC

Running from CLI with fasta file:

```bash
evaluate_dataset \
  --input example_datasets/G4_positives.fasta example_datasets/G4_negatives.fasta \
  --format fasta \
  --out_folder example_outputs/G4_dataset
```

Running from Python with CSV file with multiple sequence columns:

```python
from genbenchQC import evaluate_dataset

evaluate_dataset.run(
  inputs=['example_datasets/miRNA_mRNA_pairs_dataset.tsv'], 
  input_format='tsv', 
  out_folder='example_outputs/miRNA_mRNA_dataset', 
  sequence_column=['gene', 'noncodingRNA'], 
  label_column='label', 
  label_list=['0', '1']
)
```

**Outputs**

Running the above commands will create an output folder `example_outputs` containing the results of the evaluation. 
The output will contain (in the brackets, you can find the specific names when executing with the G4 example dataset):
- *dataset_report_\*.csv (dataset_report_label_G4_positives_vs_G4_negatives.csv)* - a CSV file with the results (pass/fail) of the evaluation between two classes of the dataset. 
- *dataset_report_\*.html (dataset_report_label_G4_positives_vs_G4_negatives.html)* - a report with the results of the evaluation between two classes of the dataset.
It contains basic information about the dataset and plots visualizing individual features of the dataset.
Each plot can be toned with red color, meaning the specific feature was too different between the two dataset classes and should be examined.
- *\*_plots (dataset_report_label_G4_positives_vs_G4_negatives_plots)* - folder containing all plots from the report in *.png* format
- *dataset_report_\*_duplicates.txt* - present only if the dataset contains duplicate sequences between classes. It contains a list of all the duplicate sequences.

For the input `example_datasets/miRNA_mRNA_pairs_dataset.tsv` where two sequence collumns are provided (`gene` and `noncodingRNA`), all report files will be generate 3 times - for sequences in`gene` column, for sequences in `noncodingRNA` column and for sequences created by concatenating `gene` and `noncodingRNA` columns.

### Split QC

```bash
evaluate_split \
  --train_input example_datasets/enhancers_train.csv \
  --test_input example_datasets/enhancers_test.csv \
  --format csv \
  --sequence_column sequence \
  --out_folder example_outputs/enhancers_dataset
```

The same evaluation would be executed from Python as follows:

```python
from genbenchQC import evaluate_split

evaluate_split.run(
  train_files=['example_datasets/enhancers_train.csv'],
  test_files=['example_datasets/enhancers_test.csv'],
  input_format='csv',
  sequence_column='sequence',
  out_folder='example_outputs/enhancers_dataset'
)
```

**Outputs**
- *train_test_check_\*_report.html (train-test_check_enhancers_train_vs_enhancers_test_report.html)* - report with clusters containing sequences from both train and test dataset parts, that could indicate data leakage.
- *train_test_check_\*.csv (train-test_check_enhancers_train_vs_enhancers_test_.csv)* - a simple report if data leakage check passed or failed.

## Input file formats

You can choose to run the tool while having different dataset formats:
- **FASTA**: The input is a FASTA file / list of FASTA files. One file needs to contain sequences of one class if running *evaluate_sequences* mode.
- **CSV/TSV**: The input is a CSV/TSV file, and you provide the name of the column containing sequences. You can have either:
  - **multiple files**, each one containing sequences from one class (similar as with FASTA input)
  - **one file** containing sequences from multiple classes. In this case, when running *evaluate_sequences* mode, you need to provide the name of the column containing class labels so the tool can split the dataset into parts. The label classes can then be inferred, or you can specify their list by yourself. The dataset will then be split into pieces containing sequences with corresponding labels and analysis will be performed similarly as with multiple files.

When having CSV/TSV input, you can also decide to provide multiple sequence columns to analyze. In this case, the analysis in modes *evaluate_sequences* and *evaluate_dataset* will be performed for each column separately and lastly for sequences made by concatenating sequences throughout all the columns. 
*evaluate_split* mode will run only the concatenated sequences.


## Development

If you want to help with the development of Genomic Benchmarks QC, you are more than welcome to join in!

For a guidance, have a look at [CONTRIBUTING.md](CONTRIBUTING.md)

## License

Genomic Benchmarks QC is MIT-style licensed, as found in the LICENSE file.