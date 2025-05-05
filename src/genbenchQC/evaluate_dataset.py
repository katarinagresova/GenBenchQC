import argparse
import pandas as pd
from pathlib import Path

from genbenchQC.utils.statistics import SequenceStatistics
from genbenchQC.utils.testing import flag_significant_differences
from genbenchQC.report.report_generator import generate_text_report, generate_html_report, generate_simple_report, generate_dataset_html_report
from genbenchQC.utils.fasta_utils import read_fasta, read_sequences_from_df, read_multisequence_df


def run_analysis(input_statistics, out_folder):
    # TODO proper naming of output files, probably accepts a list of sequences

    seqStats = SequenceStatistics(positive_fasta)
    positive_stats = seqStats.compute()
    seqStats = SequenceStatistics(negative_fasta)
    negative_stats = seqStats.compute()
    
    results = flag_significant_differences(positive_fasta, positive_stats, negative_fasta, negative_stats)

    # TODO report generation with filename, label, seq col
    txt_report_positive_path = out_folder / Path(positive_fasta.stem + '_report.txt')
    txt_report_negative_path = out_folder / Path(negative_fasta.stem + '_report.txt')
    html_report_positive_path = out_folder / Path(positive_fasta.stem + '_report.html')
    html_report_negative_path = out_folder / Path(negative_fasta.stem + '_report.html')
    simple_report_path = out_folder / 'dataset_report_simple.txt'
    html_report_path = out_folder / 'dataset_report.html'
    
    generate_text_report(positive_stats, txt_report_positive_path)
    generate_text_report(negative_stats, txt_report_negative_path)
    generate_html_report(positive_stats, html_report_positive_path)
    generate_html_report(negative_stats, html_report_negative_path)
    generate_simple_report(results, simple_report_path)
    generate_dataset_html_report(positive_stats, negative_stats, results, html_report_path)

def parse_args():
    parser = argparse.ArgumentParser(description='Evaluate positive and negative sequence datasets.')
    parser.add_argument('--input', type=str, help='Path to the dataset file. '
                                                  'Can be a list of files, each containing sequences from one class.', nargs='+', required=True)
    parser.add_argument('--format', help="Format of the input files.", choices=['fasta', 'csv', 'tsv'], default='fasta') # potentially add HF support
    parser.add_argument('--sequence_column', type=str, help='Name of the columns with sequences to analyze for datasets in CSV/TSV format. '
                                                            'Either one column or list of columns.', nargs='+', default=['sequence'])
    parser.add_argument('--label_column', type=str, help='Name with the label column for datasets in CSV/TSV format.', default='label')
    parser.add_argument('--label_list', type=str, nargs='+', help='List of label classes to consider or "infer" to parse different labels automatically from label column.'
                                                       ' For datasets in CSV/TSV format.', default=['infer'])
    parser.add_argument('--out_folder', type=str, help='Path to the output folder.', default='.')
    args = parser.parse_args()

    if args.format == 'fasta' and len(args.input) < 2:
        parser.error("When format is 'fasta', the input must contain individual files for each class.")

    return args

def main():
    try:
        args = parse_args()
    except Exception as e:
        print(f"Error parsing arguments: {e}")
        return

    if len(args.input) == 1 and args.format != 'fasta':
        # parse the input file
        if args.format == 'tsv':
            delim = '\t'
        else:
            delim = ','
        df = pd.read_csv(args.input[0], delimiter=delim)

        if args.label_column not in df.columns:
            raise ValueError(f"Label column '{args.label_column}' not found in the input file '{args.input}'.")

        # get the list of labels to consider
        if len(args.label_list) == 1 and args.label_list[0] == 'infer':
            labels = df[args.label_column].unique().tolist()
        else:
            labels = args.label_list

        for seq_col in args.sequence_column:
            seq_stats = []
            for label in labels:
                try:
                    sequences = read_sequences_from_df(df, seq_col, label)
                except Exception as e:
                    print(f"Error reading sequences': {e}")
                    return
                seq_stats += SequenceStatistics(sequences, filename=args.input[0], label=label, seq_column=seq_col)
            # TODO run statistics
        # handle multiple sequence columns by concatenating sequences
        if len(args.sequence_column) > 1:
            seq_stats = []
            for label in labels:
                sequences = read_multisequence_df(df, args.sequence_column, label)
                seq_stats += SequenceStatistics(sequences, filename=args.input[0], label=label, seq_column='_'.join(args.sequence_column))
            # TODO run statistics

    else:
        # TODO loop over input files and run statistics for each pair
        pass

if __name__ == '__main__':
    main()