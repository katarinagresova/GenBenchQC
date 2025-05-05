import argparse
from pathlib import Path

from genbenchQC.utils.statistics import SequenceStatistics
from genbenchQC.utils.testing import flag_significant_differences
from genbenchQC.report.report_generator import generate_text_report, generate_html_report, generate_simple_report, generate_dataset_html_report


def run_analysis(positive_fasta, negative_fasta, out_folder):
    # TODO proper naming of output files, probably accepts a list of sequence "files" and a list of corresponding labels

    seqStats = SequenceStatistics(positive_fasta)
    positive_stats = seqStats.compute()
    seqStats = SequenceStatistics(negative_fasta)
    negative_stats = seqStats.compute()
    
    results = flag_significant_differences(positive_fasta, positive_stats, negative_fasta, negative_stats)

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
                                                  'Can be a coma-separated list of files, each containing sequences from one class.', nargs='*', required=True)
    parser.add_argument('--format', help="Format of the input files.", choices=['fasta', 'csv', 'tsv'], default='fasta') # potentially add HF support
    parser.add_argument('--sequence_column', type=str, help='Name of the columns with sequences to analyze for datasets in CSV/TSV format. '
                                                            'Either one column or coma-separated list of columns.', nargs='*', default=['sequence'])
    parser.add_argument('--label_column', type=str, help='Name with the label column for datasets in CSV/TSV format.', default='label')
    parser.add_argument('--label_list', type=str, nargs='*', help='Comma-separated list of label classes to consider or "infer" to parse different labels automatically from label column.'
                                                       ' For datasets in CSV/TSV format.', default=['infer'])
    parser.add_argument('--out_folder', type=str, help='Path to the output folder.', default='.')
    return parser.parse_args()

def main():
    args = parse_args()

    if len(args.input) == 1:
        if len(args.label_list) == 1 and args.label_list[0] == 'infer':
            # TODO infer labels from the input file
            pass
        else:
            # TODO check if the labels are in the input file
            pass
        # TODO loop through individual classes and run statistics for each pair
    else:
        # TODO loop over input files and run statistics for each pair
        pass

    #run_analysis(
    #    Path(args.positive_fasta),
    #    Path(args.negative_fasta),
    #    Path(args.out_folder)
    #)

if __name__ == '__main__':
    main()