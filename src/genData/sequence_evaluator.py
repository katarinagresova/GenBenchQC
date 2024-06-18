import argparse

from utils.statistics import compute_sequence_statistics, evaluate_stats
from utils.report_generator import generate_text_report, generate_html_report, generate_simple_report

def evaluate_positive_sequences(fasta_file, txt_report_path, html_report_path, simple_report_path, threshold):
    
    stats = compute_sequence_statistics(fasta_file)
    results = evaluate_stats(stats, threshold)
    
    if txt_report_path is not None:
        generate_text_report(stats, results, txt_report_path)
    if html_report_path is not None:
        generate_html_report(stats, results, html_report_path)
    if simple_report_path is not None:
        generate_simple_report(results, simple_report_path)

def parse_args():
    parser = argparse.ArgumentParser(description='Evaluate a positive sequence dataset.')
    parser.add_argument('input', type=str, help='Path to the positive dataset in FASTA format.')
    parser.add_argument('--txt_report_path', type=str, help='Path to the output text report.', default=None)
    parser.add_argument('--html_report_path', type=str, help='Path to the output HTML report.', default=None)
    parser.add_argument('--simple_report_path', type=str, help='Path to the output simple report.', default=None)
    parser.add_argument('--threshold', type=float, help='P-value threshold for flagging significant differences.', default=0.01)
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    evaluate_positive_sequences(
        args.input, 
        args.txt_report_path, 
        args.html_report_path, 
        args.simple_report_path, 
        args.threshold
    )