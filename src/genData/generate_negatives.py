import argparse
from pathlib import Path

def run(pos_fasta, neg_fasta):
    print("Not implemented yet")

def parse_args():
    parser = argparse.ArgumentParser(description='Generate negative sequences.')
    parser.add_argument('positive_fasta', type=str, help='Path to the positive dataset in FASTA format.')
    parser.add_argument('negative_fasta', type=str, help='Path to output FASTA file for generated negative sequences.')
    return parser.parse_args()

def main():
    args = parse_args()
    run(Path(args.positive_fasta), Path(args.negative_fasta))

if __name__ == '__main__':
    main()