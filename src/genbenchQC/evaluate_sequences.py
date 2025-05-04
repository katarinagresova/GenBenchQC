import argparse
from pathlib import Path

from genbenchQC.utils.statistics import SequenceStatistics
from genbenchQC.report.report_generator import generate_text_report, generate_html_report

def run(fasta_file, out_folder):
    
    seqStats = SequenceStatistics(fasta_file)
    stats = seqStats.compute()
    
    txt_report_path = Path(out_folder, str(fasta_file.stem) + '_report.txt')
    html_report_path = Path(out_folder, str(fasta_file.stem) + '_report.html')

    generate_text_report(stats, txt_report_path)
    generate_html_report(stats, html_report_path)

def parse_args():
    parser = argparse.ArgumentParser(description='A tools for evaluating sequence data.')
    parser.add_argument('input', type=str, help='Path to the fasta file.')
    parser.add_argument('--out_folder', type=str, help='Path to the output folder.', default='.')
    return parser.parse_args()

def main():
    args = parse_args()
    run(
        Path(args.input), 
        Path(args.out_folder), 
    )

if __name__ == '__main__':
    main()