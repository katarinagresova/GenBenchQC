from Bio import SeqIO, SeqRecord, Seq

def read_fasta(fasta_file):
    return [str(record.seq).upper() for record in SeqIO.parse(fasta_file, 'fasta')]

def write_fasta(sequences, output_file):
    records = [SeqRecord.SeqRecord(Seq.Seq(seq), id=f'seq_{i}') for i, seq in enumerate(sequences)]
    SeqIO.write(records, output_file, 'fasta')