# Detecting and Mitigating Biases in Genomic Data

## Installation

Clone this repository

```bash
git clone https://github.com/katarinagresova/MitigatingBiases.git
```

Switch to the directory

```bash
cd MitigatingBiases
```

Switch to `refactor` branch

```bash
git checkout refactor
```

Prepare conda environment

```bash
. prepare_conda.sh
```


## Example

Prepare sample data

```bash
python sample_data/generate_fasta.py
```

Run QC tool on positive data

```bash
python src/genData/sequence_evaluator sample_data/pos.fasta --txt_report_path sample_data/pos_report.txt --simple_report_path sample_data/pos_simple_report.txt
```
