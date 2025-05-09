from Bio import SeqIO, SeqRecord, Seq

def read_fasta(fasta_file):
    return [str(record.seq).upper() for record in SeqIO.parse(fasta_file, 'fasta')]

def write_fasta(sequences, output_file):
    records = [SeqRecord.SeqRecord(Seq.Seq(seq), id=f'seq_{i}') for i, seq in enumerate(sequences)]
    SeqIO.write(records, output_file, 'fasta')

def read_sequences_from_df(df, seq_column, label_column=None, label=None):
    if seq_column not in df.columns:
        raise ValueError(f"Sequence column '{seq_column}' not found in the input file.")

    if label_column is None:
        return [str(seq).upper() for seq in df[seq_column].tolist()]

    df_parsed = df[df[label_column] == label]

    if df_parsed.empty:
        raise ValueError(f"No sequences found for label '{label}' in column '{label_column}'.")

    return [str(seq).upper() for seq in df_parsed[seq_column].tolist()]

def read_multisequence_df(df, seq_columns, label_column=None, label=None):
    all_sequences = [read_sequences_from_df(df, seq_column, label_column, label) for seq_column in seq_columns]
    concatenated_sequences = [''.join(seqs) for seqs in zip(*all_sequences)]
    return concatenated_sequences


