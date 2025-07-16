import argparse
import logging
from pathlib import Path
import os
from cdhit_reader import read_cdhit

from genbenchQC.utils.input_utils import setup_logger, read_files_to_sequence_list, write_fasta

def run_clustering(fasta_file, clustered_file):
    logging.info("Running CD-HIT clustering.")
    os.system(f"cd-hit-est -i {fasta_file} -o {clustered_file} -c 0.95 -n 10 -aS 0.8 -aL 0.8")
    clusters = read_cdhit(f"{clustered_file}.clstr").read_items()
    logging.debug(f"CD-HIT clustering completed. {len(clusters)} clusters found.")

    # Collect clusters with >1 sequence and mixed train/test entries
    mixed_clusters = []

    for cluster in clusters:
        if len(cluster.sequences) == 1:
            continue
        seq_ids = [seq.name for seq in cluster.sequences]
        has_train = any("_train" in seq_id for seq_id in seq_ids)
        has_test = any("_test" in seq_id for seq_id in seq_ids)

        if has_train and has_test:
            mixed_clusters.append(seq_ids)
    return mixed_clusters

def run(train_files, test_files, input_format, out_folder, sequence_column):

    logging.info("Starting train-test split evaluation.")

    if not Path(out_folder).exists():
        logging.info(f"Output folder {out_folder} does not exist. Creating it.")
        Path(out_folder).mkdir(parents=True, exist_ok=True)

    Path(out_folder, "tmp").mkdir(parents=True, exist_ok=True)

    train_sequences = read_files_to_sequence_list(train_files, input_format, sequence_column)
    train_index = [f"{i}_train" for i in range(len(train_sequences))]
    logging.info(f"Read {len(train_sequences)} sequences from training files.")

    test_sequences = read_files_to_sequence_list(test_files, input_format, sequence_column)
    test_index = [f"{i}_test" for i in range(len(test_sequences))]
    logging.info(f"Read {len(test_sequences)} sequences from testing files.")

    sequences = train_sequences + test_sequences
    index = train_index + test_index
    merged_fasta_path = Path(out_folder, "tmp") / 'train_test_sequences.fasta'
    write_fasta(sequences, merged_fasta_path, index)

    clusters = run_clustering(merged_fasta_path, Path(out_folder, "tmp/clustered_sequences"))
    logging.debug(f"Having {len(clusters)} mixed clusters: {clusters}")

    # TODO delete tmp folder

    logging.info("Train-test split evaluation successfully completed.")

def parse_args():
    parser = argparse.ArgumentParser(description='Check data leakage in dataset train-test split.')
    parser.add_argument('--train_input', type=str, help='Path to the dataset file with training data. Can be multiple files that will be evaluated as one dataset part.', nargs='+', required=True)
    parser.add_argument('--test_input', type=str, help='Path to the dataset file with testing data. Can be multiple files that will be evaluated as one dataset part.', nargs='+',
                        required=True)
    parser.add_argument('--format', help="Format of the input files.", choices=['fasta', 'csv', 'tsv'], required=True)
    parser.add_argument('--sequence_column', type=str, help='Name of the columns with sequences to analyze for datasets in CSV/TSV format. '
                                                            'Either one column or list of columns.', nargs='+', default=['sequence'])
    parser.add_argument('--out_folder', type=str, help='Path to the output folder.', default='.')
    parser.add_argument('--log_level', type=str, help='Logging level, default to INFO.', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='INFO')
    parser.add_argument('--log_file', type=str, help='Path to the log file.', default=None)
    args = parser.parse_args()

    return args

def main():
    args = parse_args()
    setup_logger(args.log_level, args.log_file)
    run(args.train_input, args.test_input, args.format, args.out_folder, args.sequence_column)

if __name__ == '__main__':
    main()