import argparse
from pathlib import Path
from typing import Optional
import logging

from genbenchQC.utils.statistics import SequenceStatistics
from genbenchQC.report.report_generator import generate_json_report, generate_sequence_html_report
from genbenchQC.utils.input_utils import read_fasta, read_sequences_from_df, read_multisequence_df, read_csv_file, setup_logger

def run_analysis(seq_stats, out_folder, report_types):

    if not Path(out_folder).exists():
        logging.info(f"Output folder {out_folder} does not exist. Creating it.")
        Path(out_folder).mkdir(parents=True, exist_ok=True)

    stats, end_position = seq_stats.compute()

    filename = Path(seq_stats.filename).stem
    if seq_stats.seq_column is not None:
        filename += f'_{seq_stats.seq_column}'
    if seq_stats.label is not None:
        filename += f'_{seq_stats.label}'

    if 'json' in report_types:
        json_report_path = Path(out_folder, filename + '_report.json')
        generate_json_report(stats, json_report_path)
    if 'html' in report_types:
        html_report_path = Path(out_folder, filename + '_report.html')
        plots_path = out_folder / Path(filename + '_plots')
        generate_sequence_html_report(stats, html_report_path, plots_path, end_position=end_position)

def run(input_file, 
        input_format, 
        out_folder='.', 
        sequence_column: Optional[list[str]] = ['sequences'], 
        label_column=None, 
        label: Optional[str] = None,
        report_types: Optional[list[str]] = ['html'],
        end_position: Optional[int] = None):
    
    logging.info("Starting sequence evaluation.")

    if input_format == 'fasta':
        seqs = read_fasta(input_file)
        logging.debug(f"Read {len(seqs)} sequences from FASTA file.")
        run_analysis(
            SequenceStatistics(seqs, Path(input_file).name, label=label),
            out_folder, report_types=report_types, end_position=end_position
        )
    else:
        df = read_csv_file(input_file, input_format, sequence_column, label_column)

        for seq_col in sequence_column:
            sequences = read_sequences_from_df(df, seq_col, label_column, label)
            logging.debug(f"Read {len(sequences)} sequences from CSV/TSV file.")
            run_analysis(
                SequenceStatistics(sequences, filename=Path(input_file).name, seq_column=seq_col, label=label), 
                out_folder, report_types=report_types, end_position=end_position
            )

        if len(sequence_column) > 1:
            sequences = read_multisequence_df(df, sequence_column, label_column, label)
            run_analysis(
                SequenceStatistics(sequences, filename=Path(input_file).name, seq_column='_'.join(sequence_column), label=label), 
                out_folder, report_types=report_types, end_position=end_position
            )

    logging.info("Sequence evaluation successfully completed.")

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
    parser.add_argument('--report_types', type=str, nargs='+', choices=['json', 'html'],
                        help='Types of reports to generate. Default: [html]', default=['html'])
    parser.add_argument('--end_position', type=int, help='End position of the sequences to plot in the per position plots.', default=None)
    parser.add_argument('--log_level', type=str, help='Logging level, default to INFO.',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='INFO')
    parser.add_argument('--log_file', type=str, help='Path to the log file.', default=None)

    args = parser.parse_args()

    if (args.label_column is not None and args.label is None) or (args.label_column is None and args.label is not None):
        parser.error("--label_column and --label must be provided together.")

    return args

def main():
    args = parse_args()
    setup_logger(args.log_level, args.log_file)
    run(args.input, 
        args.format, 
        args.out_folder, 
        args.sequence_column, 
        args.label_column, 
        args.label, 
        args.report_types,
        args.end_position
    )

if __name__ == '__main__':
    main()