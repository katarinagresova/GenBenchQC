import argparse
from pathlib import Path
from typing import Optional

from genbenchQC.utils.statistics import SequenceStatistics
from genbenchQC.report.report_generator import generate_text_report, generate_html_report
from genbenchQC.utils.input_utils import read_fasta, read_sequences_from_df, read_multisequence_df, read_csv_file

def run_analysis(seq_stats, out_folder):
    stats = seq_stats.compute()

    filename = Path(seq_stats.filename).stem
    if seq_stats.seq_column is not None:
        filename += f'_{seq_stats.seq_column}'

    txt_report_path = Path(out_folder, filename + '_report.txt')
    html_report_path = Path(out_folder, filename + '_report.html')

    generate_text_report(stats, txt_report_path)
    generate_html_report(stats, html_report_path)

def run(input_file, input_format, out_folder='.', sequence_column: Optional[list[str]] = ['sequences']):
    if input_format == 'fasta':
        seqs = read_fasta(input_file)
        run_analysis(SequenceStatistics(seqs, input_file), out_folder)
    else:
        df = read_csv_file(input_file, input_format, sequence_column)

        for seq_col in sequence_column:
            sequences = read_sequences_from_df(df, seq_col)
            run_analysis(SequenceStatistics(sequences, filename=input_file, seq_column=seq_col), out_folder)

        if len(sequence_column) > 1:
            sequences = read_multisequence_df(df, sequence_column)
            run_analysis(SequenceStatistics(sequences, filename=input_file, seq_column='_'.join(sequence_column)), out_folder)

def parse_args():
    parser = argparse.ArgumentParser(description='A tools for evaluating sequence data.')
    parser.add_argument('--input', type=str, help='Path to the input file.', required=True)
    parser.add_argument('--format', help="Format of the input file.", choices=['fasta', 'csv', 'tsv'], required=True)
    parser.add_argument('--sequence_column', type=str,
                        help='Name of the columns with sequences to analyze for datasets in CSV/TSV format. '
                             'Either one column or list of columns.', nargs='+', default=['sequence'])
    parser.add_argument('--out_folder', type=str, help='Path to the output folder.', default='.')
    return parser.parse_args()

def main():
    args = parse_args()
    run(args.input, args.format, args.out_folder, args.sequence_column)


if __name__ == '__main__':
    main()