from Bio import SeqIO, SeqRecord, Seq
import pandas as pd

def read_fasta(fasta_file):
    return [str(record.seq).upper() for record in SeqIO.parse(fasta_file, 'fasta')]

def write_fasta(sequences, output_file):
    records = [SeqRecord.SeqRecord(Seq.Seq(seq), id=f'seq_{i}') for i, seq in enumerate(sequences)]
    SeqIO.write(records, output_file, 'fasta')

def read_csv_file(file_path, input_format, seq_columns, label_columns=None):
    delim = '\t' if input_format == 'tsv' else ','

    columns = seq_columns.copy()
    if label_columns is not None:
        columns += [label_columns]
    df = pd.read_csv(file_path, delimiter=delim, usecols=columns)
    df[seq_columns] = df[seq_columns].apply(lambda col: col.str.upper())

    return df


def read_sequences_from_df(df, seq_column, label_column=None, label=None):
    if label_column is None:
        return df[seq_column].tolist()

    df_parsed = df[df[label_column] == label]

    if df_parsed.empty:
        raise ValueError(f"No sequences found for label '{label}' in column '{label_column}'.")

    return df_parsed[seq_column].tolist()

def read_multisequence_df(df, seq_columns, label_column=None, label=None):
    all_sequences = [read_sequences_from_df(df, seq_column, label_column, label) for seq_column in seq_columns]
    concatenated_sequences = [''.join(seqs) for seqs in zip(*all_sequences)]
    return concatenated_sequences


