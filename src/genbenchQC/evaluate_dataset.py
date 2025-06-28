import argparse
import logging
from pathlib import Path
from itertools import combinations
from typing import Optional

from genbenchQC.utils.statistics import SequenceStatistics
from genbenchQC.utils.testing import flag_significant_differences
from genbenchQC.report.report_generator import generate_json_report, generate_sequence_html_report, generate_simple_report, generate_dataset_html_report
from genbenchQC.utils.input_utils import read_fasta, read_sequences_from_df, read_multisequence_df, read_csv_file, setup_logger

def run_analysis(input_statistics, out_folder, report_types, seq_report_types, threshold=0.015):
   
    out_folder = Path(out_folder)

    # run individual analysis
    for s in input_statistics:
        stats, end_position = s.compute()

        if seq_report_types:

            filename = Path(s.filename).stem
            if s.seq_column is not None:
                filename += f'_{s.seq_column}'
            if s.label is not None:
                filename += f'_{s.label}'

            if 'json' in seq_report_types:
                json_report_path = out_folder / Path(filename + '_report.json')
                generate_json_report(stats, json_report_path)

            if 'html' in seq_report_types:
                html_report_path = out_folder / Path(filename + '_report.html')
                plots_path = out_folder / Path(filename + '_plots')
                generate_sequence_html_report(stats, html_report_path, plots_path, end_position)

    if len(input_statistics) < 2:
        return

    # run pair comparison analysis with all combinations
    for stat1, stat2 in combinations(input_statistics, 2):
        filename = "dataset_report"
        if stat1.seq_column is not None:
            filename += f'_{stat1.seq_column}'
        if stat1.label is not None and stat2.label is not None:
            filename += f'_label_{stat1.label}_vs_{stat2.label}'
            logging.debug(f"Comparing datasets label: {stat1.label} vs {stat2.label}")
        else:
            filename += f'_{Path(stat1.filename).stem}_{Path(stat2.filename).stem}'
            logging.debug(
                f"Comparing datasets: {stat1.filename} vs {stat2.filename}")

        results = flag_significant_differences(
            stat1.sequences, stat1.stats, 
            stat2.sequences, stat2.stats, 
            threshold=threshold, 
            end_position=min(stat1.end_position, stat2.end_position)
        )
        
        if 'simple' in report_types:
            simple_report_path = out_folder / Path(f'{filename}.csv')
            generate_simple_report(results, simple_report_path)

        if 'html' in report_types:
            html_report_path = out_folder / Path(f'{filename}.html')
            plots_path = out_folder / Path(f'{filename}_plots')
            generate_dataset_html_report(
                stat1, stat2, results, 
                html_report_path, 
                plots_path=plots_path, 
                threshold=threshold,
                end_position=min(stat1.end_position, stat2.end_position)
            )

def run(inputs, 
        input_format, 
        out_folder='.', 
        sequence_column: Optional[list[str]] = ['sequences'], 
        label_column='label', 
        label_list: Optional[list[str]] = ['infer'],
        regression: Optional[bool] = False,
        report_types: Optional[list[str]] = ['html', 'simple'],
        seq_report_types: Optional[list[str]] = None,
        end_position: Optional[int] = None):

    logging.info("Starting dataset evaluation.")

    if not Path(out_folder).exists():
        logging.info(f"Output folder {out_folder} does not exist. Creating it.")
        Path(out_folder).mkdir(parents=True, exist_ok=True)

    # we have multiple fasta files with one label each
    if input_format == 'fasta':
        seq_stats = []
        for input_file in inputs:
            sequences = read_fasta(input_file)
            logging.debug(f"Read {len(sequences)} sequences from FASTA file {input_file}.")
            seq_stats += [SequenceStatistics(sequences, filename=Path(input_file).name, 
                                             label=Path(input_file).stem, end_position=end_position)]
        run_analysis(seq_stats, out_folder, report_types, seq_report_types)

    # we have CSV/TSV
    else:
        # we have one file with multiple labels or regression target
        if len(inputs) == 1:
            df = read_csv_file(inputs[0], input_format, sequence_column, label_column)

            # if regression is True, we split the label column into two classes
            if regression:
                # infer the threshold as the median of the label column
                threshold = df[label_column].median()
                logging.debug(f"Inferred threshold for regression: {threshold}")
                df[label_column] = df[label_column].apply(lambda x: 'high' if x >= threshold else 'low')
                labels = ['high', 'low']

            # get the list of labels to consider
            elif len(label_list) == 1 and label_list[0] == 'infer':
                labels = df[label_column].unique().tolist()
                logging.debug(f"Inferred labels: {labels}")
            else:
                labels = label_list

            # loop over sequences with specific label and run statistics
            for seq_col in sequence_column:
                seq_stats = []
                for label in labels:
                    sequences = read_sequences_from_df(df, seq_col, label_column, label)
                    logging.debug(f"Read {len(sequences)} sequences for label '{label}' from column '{seq_col}'.")
                    seq_stats += [SequenceStatistics(sequences, filename=Path(inputs[0]).name, label=label, 
                                                     seq_column=seq_col, end_position=end_position)]
                run_analysis(seq_stats, out_folder, report_types, seq_report_types)

            # handle multiple sequence columns by concatenating sequences and running statistics on them
            if len(sequence_column) > 1:
                seq_stats = []
                for label in labels:
                    sequences = read_multisequence_df(df, sequence_column, label_column, label)
                    seq_stats += [SequenceStatistics(sequences, filename=Path(inputs[0]).name, label=label,
                                                     seq_column='_'.join(sequence_column))]
                run_analysis(seq_stats, out_folder, report_types, seq_report_types)

        # we have multiple files with one label each
        else:
            # run statistics across input files
            for seq_col in sequence_column:
                seq_stats = []
                for input_file in inputs:
                    sequences = read_sequences_from_df(read_csv_file(input_file, input_format, seq_col), seq_col)
                    logging.debug(f"Read {len(sequences)} sequences from file {input_file} in column '{seq_col}'.")
                    seq_stats += [SequenceStatistics(sequences, filename=Path(input_file).name, 
                                                     label=Path(input_file).stem, seq_column=seq_col,
                                                     end_position=end_position)]
                run_analysis(seq_stats, out_folder, report_types, seq_report_types)

            # handle multiple sequence columns
            if len(sequence_column) > 1:
                seq_stats = []
                for input_file in inputs:
                    sequences = read_multisequence_df(read_csv_file(input_file, input_format, sequence_column), sequence_column)
                    seq_stats += [SequenceStatistics(sequences, filename=Path(input_file).name, label=Path(input_file).stem,
                                                     seq_column='_'.join(sequence_column), end_position=end_position)]
                run_analysis(seq_stats, out_folder, report_types, seq_report_types)

    logging.info("Dataset evaluation successfully completed.")


def parse_args():
    parser = argparse.ArgumentParser(description='Evaluate positive and negative sequence datasets.')
    parser.add_argument('--input', type=str, help='Path to the dataset file. '
                                                  'Can be a list of files, each containing sequences from one class.', nargs='+', required=True)
    parser.add_argument('--format', help="Format of the input files.", choices=['fasta', 'csv', 'tsv'], required=True) # potentially add HF support
    parser.add_argument('--sequence_column', type=str, help='Name of the columns with sequences to analyze for datasets in CSV/TSV format. '
                                                            'Either one column or list of columns.', nargs='+', default=['sequence'])
    parser.add_argument('--label_column', type=str, help='Name with the label column for datasets in CSV/TSV format.', default='label')
    parser.add_argument('--label_list', type=str, nargs='+', help='List of label classes to consider or "infer" to parse different labels automatically from label column.'
                                                       ' For datasets in CSV/TSV format.', default=['infer'])
    parser.add_argument('--regression', action='store_true', help='If True, label column is considered as a regression target and values are split into 2 classes')
    parser.add_argument('--out_folder', type=str, help='Path to the output folder.', default='.')
    parser.add_argument('--report_types', type=str, nargs='+', choices=['json', 'html', 'simple'], default=['html', 'simple'],
                        help='Types of reports to generate. Default: [html, simple].')
    parser.add_argument('--seq_report_types', type=str, nargs='+', choices=['json', 'html'], default=[],
                        help='Types of reports to generate for individual groups of sequences. Default: [].')
    parser.add_argument('--end_position', type=int, help='End position of the sequences to consider in per position statistics.', default=None)
    parser.add_argument('--log_level', type=str, help='Logging level, default to INFO.', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='INFO')
    parser.add_argument('--log_file', type=str, help='Path to the log file.', default=None)
    args = parser.parse_args()

    if args.format == 'fasta' and len(args.input) < 2:
        parser.error("When format is 'fasta', the input must contain individual files for each class.")

    return args

def main():
    args = parse_args()
    setup_logger(args.log_level, args.log_file)
    run(args.input, 
        args.format, 
        args.out_folder, 
        args.sequence_column, 
        args.label_column, 
        args.label_list,
        args.regression,
        args.report_types,
        args.seq_report_types,
        args.end_position
    )

if __name__ == '__main__':
    main()