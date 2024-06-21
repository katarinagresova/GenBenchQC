import argparse
from pathlib import Path

from genData.utils.statistics import compute_sequence_statistics
from genData.utils.testing import flag_significant_differences
from genData.utils.report_generator import generate_text_report, generate_html_report, generate_simple_report, generate_dataset_html_report

def run(positive_fasta, negative_fasta, out_folder):
    
    positive_stats = compute_sequence_statistics(positive_fasta)
    negative_stats = compute_sequence_statistics(negative_fasta)
    
    results = flag_significant_differences(positive_fasta, positive_stats, negative_fasta, negative_stats)

    txt_report_positive_path = out_folder / Path(positive_fasta.stem + '_report.txt')
    txt_report_negative_path = out_folder / Path(negative_fasta.stem + '_report.txt')
    html_report_positive_path = out_folder / Path(positive_fasta.stem + '_report.html')
    html_report_negative_path = out_folder / Path(negative_fasta.stem + '_report.html')
    simple_report_path = out_folder / 'dataset_report_simple.txt'
    html_report_path = out_folder / 'dataset_report.html'
    
    generate_text_report(positive_stats, txt_report_positive_path)
    generate_text_report(negative_stats, txt_report_negative_path)
    generate_html_report(positive_stats, results, html_report_positive_path)
    generate_html_report(negative_stats, results, html_report_negative_path)
    generate_simple_report(results, simple_report_path)
    generate_dataset_html_report(positive_stats, negative_stats, results, html_report_path)

def parse_args():
    parser = argparse.ArgumentParser(description='Evaluate positive and negative sequence datasets.')
    parser.add_argument('positive_fasta', type=str, help='Path to the positive dataset in FASTA format.')
    parser.add_argument('negative_fasta', type=str, help='Path to the negative dataset in FASTA format.')
    parser.add_argument('--out_folder', type=str, help='Path to the output folder.', default='.')
    return parser.parse_args()

def main():
    args = parse_args()
    run(
        Path(args.positive_fasta), 
        Path(args.negative_fasta), 
        Path(args.out_folder)
    )

if __name__ == '__main__':
    main()