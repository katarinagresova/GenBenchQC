from Bio import SeqIO, SeqRecord, Seq
import pandas as pd
import logging
import json

def read_fasta(fasta_file):
    logging.debug(f"Reading FASTA file: {fasta_file}")
    return [str(record.seq).upper() for record in SeqIO.parse(fasta_file, 'fasta')]

def write_fasta(sequences, output_file):
    records = [SeqRecord.SeqRecord(Seq.Seq(seq), id=f'seq_{i}') for i, seq in enumerate(sequences)]
    SeqIO.write(records, output_file, 'fasta')

def read_csv_file(file_path, input_format, seq_columns, label_columns=None):
    delim = '\t' if input_format == 'tsv' else ','

    columns = seq_columns.copy()
    if label_columns is not None:
        columns += [label_columns]
    df = pd.read_csv(file_path, delimiter=delim, usecols=columns, dtype=str)
    df[seq_columns] = df[seq_columns].apply(lambda col: col.str.upper())

    logging.debug(f"Read CSV/TSV file: {file_path}, shape: {df.shape}, columns: {columns}")

    return df

def read_sequences_from_df(df, seq_column, label_column=None, label=None):
    if label_column is None:
        return df[seq_column].tolist()

    logging.debug(f"Filtering sequences by label: {label} in column: {label_column}")
    df_parsed = df[df[label_column] == label]

    if df_parsed.empty:
        logging.error(f"No sequences found for label '{label}' in column '{label_column}'.")
        raise ValueError(f"No sequences found for label '{label}' in column '{label_column}'.")

    return df_parsed[seq_column].tolist()

def read_multisequence_df(df, seq_columns, label_column=None, label=None):
    logging.debug(f"Concatenating sequences from multiple columns: {label_column}")
    all_sequences = [read_sequences_from_df(df, seq_column, label_column, label) for seq_column in seq_columns]
    concatenated_sequences = [''.join(seqs) for seqs in zip(*all_sequences)]
    return concatenated_sequences

def setup_logger(level=logging.INFO, file=None):
    if file:
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(file, mode='w'),
                logging.StreamHandler()
            ]
        )
    else:
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.StreamHandler()
            ]
        )

    # Suppress matplotlib debug messages
    logging.getLogger("matplotlib").setLevel(logging.WARNING)

def read_stats_json(stats_json_file):
    with open(stats_json_file, 'r') as file:
        stats = json.load(file)
    return stats

def write_stats_json(stats, stats_json_file):
    with open(stats_json_file, 'w') as file:
        json.dump(stats, file, indent=4)
