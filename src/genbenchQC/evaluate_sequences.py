import argparse
from pathlib import Path
import pandas as pd

from genbenchQC.utils.statistics import SequenceStatistics
from genbenchQC.report.report_generator import generate_text_report, generate_html_report
from genbenchQC.utils.fasta_utils import read_fasta, read_sequences_from_df, read_multisequence_df

def run(seq_stats, out_folder):
    stats = seq_stats.compute()
    
    txt_report_path = Path(out_folder, str(Path(seq_stats.filename).stem) + '_report.txt')
    html_report_path = Path(out_folder, str(Path(seq_stats.filename).stem) + '_report.html')

    generate_text_report(stats, txt_report_path)
    generate_html_report(stats, html_report_path)

def parse_args():
    parser = argparse.ArgumentParser(description='A tools for evaluating sequence data.')
    parser.add_argument('input', type=str, help='Path to the input file.', required=True)
    parser.add_argument('--format', help="Format of the input file.", choices=['fasta', 'csv', 'tsv'], required=True)
    parser.add_argument('--sequence_column', type=str,
                        help='Name of the columns with sequences to analyze for datasets in CSV/TSV format. '
                             'Either one column or list of columns.', nargs='+', default=['sequence'])
    parser.add_argument('--out_folder', type=str, help='Path to the output folder.', default='.')
    return parser.parse_args()

def main():
    args = parse_args()

    if args.format == 'fasta':
        seqs = read_fasta(args.input)
        run(SequenceStatistics(seqs, args.input), args.out_folder)
    else:
        sep = ',' if args.format == 'csv' else '\t'
        df = pd.read_csv(args.input, delimiter=sep)

        for seq_col in args.sequence_column:
            try:
                sequences = read_sequences_from_df(df, seq_col)
            except Exception as e:
                print(f"Error reading sequences: {e}")
                return
            run(SequenceStatistics(sequences, filename=args.input, seq_column=seq_col), args.out_folder)

        if len(args.sequence_column) > 1:
            sequences = read_multisequence_df(df, args.sequence_column)
            run(SequenceStatistics(sequences, filename=args.input, seq_column='_'.join(args.sequence_column)), args.out_folder)


if __name__ == '__main__':
    main()