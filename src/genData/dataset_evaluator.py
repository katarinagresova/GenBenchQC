import argparse

from utils.fasta_utils import read_fasta
from utils.statistics import compute_sequence_statistics, flag_significant_differences
from utils.report_generator import generate_text_report, generate_html_report

def evaluate_datasets(positive_fasta, negative_fasta, txt_report_path, html_report_path):
    positive_sequences = read_fasta(positive_fasta)
    negative_sequences = read_fasta(negative_fasta)
    
    positive_stats = compute_sequence_statistics(positive_sequences)
    negative_stats = compute_sequence_statistics(negative_sequences)
    
    flagged_stats = flag_significant_differences(positive_stats, negative_stats)
    
    txt_report_path = generate_text_report(flagged_stats, txt_report_path)
    html_report_path = generate_html_report(flagged_stats, html_report_path)

def parse_args():
    parser = argparse.ArgumentParser(description='Evaluate positive and negative sequence datasets.')
    parser.add_argument('positive_fasta', type=str, help='Path to the positive dataset in FASTA format.')
    parser.add_argument('negative_fasta', type=str, help='Path to the negative dataset in FASTA format.')
    parser.add_argument('txt_report_path', type=str, help='Path to the output text report.')
    parser.add_argument('html_report_path', type=str, help='Path to the output HTML report.')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    evaluate_datasets(args.positive_fasta, args.negative_fasta, args.txt_report_path, args.html_report_path)