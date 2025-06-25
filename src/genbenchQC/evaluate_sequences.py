import argparse
from pathlib import Path
from typing import Optional

from genbenchQC.utils.statistics import SequenceStatistics
from genbenchQC.report.report_generator import generate_json_report, generate_html_report
from genbenchQC.utils.input_utils import read_fasta, read_sequences_from_df, read_multisequence_df, read_csv_file

def run_analysis(seq_stats, out_folder):

    if not Path(out_folder).exists():
        print(f"Output folder {out_folder} does not exist. Creating it.")
        Path(out_folder).mkdir(parents=True, exist_ok=True)

    stats = seq_stats.compute()

    filename = Path(seq_stats.filename).stem
    if seq_stats.seq_column is not None:
        filename += f'_{seq_stats.seq_column}'
    if seq_stats.label is not None:
        filename += f'_{seq_stats.label}'

    json_report_path = Path(out_folder, filename + '_report.json')
    html_report_path = Path(out_folder, filename + '_report.html')

    generate_json_report(stats, json_report_path)
    generate_html_report(stats, html_report_path)

def run(input_file, input_format, out_folder='.', sequence_column: Optional[list[str]] = ['sequences'], label_column=None, label: Optional[str] = None):
    
    if input_format == 'fasta':
        seqs = read_fasta(input_file)
        run_analysis(
            SequenceStatistics(seqs, Path(input_file).name, label=label), 
            out_folder
        )
    else:
        df = read_csv_file(input_file, input_format, sequence_column, label_column)

        for seq_col in sequence_column:
            sequences = read_sequences_from_df(df, seq_col, label_column, label)
            run_analysis(
                SequenceStatistics(sequences, filename=Path(input_file).name, seq_column=seq_col, label=label), 
                out_folder
            )

        if len(sequence_column) > 1:
            sequences = read_multisequence_df(df, sequence_column, label_column, label)
            run_analysis(
                SequenceStatistics(sequences, filename=Path(input_file).name, seq_column='_'.join(sequence_column), label=label), 
                out_folder
            )

def parse_args():
    parser = argparse.ArgumentParser(description='A tools for evaluating sequence data.')
    parser.add_argument('--input', type=str, help='Path to the input file.', required=True)
    parser.add_argument('--format', help="Format of the input file.", choices=['fasta', 'csv', 'tsv'], required=True)
    parser.add_argument('--sequence_column', type=str,
                        help='Name of the columns with sequences to analyze for datasets in CSV/TSV format. '
                             'Either one column or list of columns.', nargs='+', default=['sequence'])
    parser.add_argument('--label_column', type=str, help='Name with the label column for datasets in CSV/TSV format. Needed only if you want to select a specific class from the dataset.',
                        default=None)
    parser.add_argument('--label', type=str,
                        help='Label of the class to select from the whole dataset. If not specified, the whole dataset is taken and analyzed as one piece.', default=None)
    parser.add_argument('--out_folder', type=str, help='Path to the output folder.', default='.')

    args = parser.parse_args()

    if (args.label_column is not None and args.label is None) or (args.label_column is None and args.label is not None):
        parser.error("--label_column and --label must be provided together.")

    return args

def main():
    args = parse_args()
    run(args.input, args.format, args.out_folder, args.sequence_column, args.label_column, args.label)


if __name__ == '__main__':
    main()