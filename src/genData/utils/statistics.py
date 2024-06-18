import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, ranksums
from statsmodels.stats.multitest import fdrcorrection
from collections import Counter

from utils.fasta_utils import read_fasta

def compute_sequence_statistics(fasta_file):
    """
    Compute various statistics from the given list of sequences.
    @param sequences: A list of sequences.
    @return: A dictionary containing the following statistics:
        - Filename: str
        - Number of sequences: int
        - Number of bases: int
        - %GC content: float
        - Per sequence nucleotide content: dict {sequence_id: dict {nucleotide: count}}
        - Per sequence dinucleotide content: dict {sequence_id: dict {dinucleotide: count}}
        - Per position nucleotide content: dict {position: dict {nucleotide: count}}
        - Per position reversed nucleotide content: dict {position: dict {nucleotide: count}}
        - Per sequence GC content: dict {sequence_id: float}
        - Sequence lenghts: dict {sequence_id: int}
        - Sequence duplication levels: dict {sequence: count}
    """

    sequences = read_fasta(fasta_file)
    stats = {
        'Filename': fasta_file,
        'Number of sequences': len(sequences),
        'Number of bases': sum(len(sequence) for sequence in sequences),
        '%GC content': sum(sequence.count('G') + sequence.count('C') for sequence in sequences) / sum(len(sequence) for sequence in sequences),
        'Per sequence nucleotide content': compute_per_sequence_nucleotide_content(sequences),
        'Per sequence dinucleotide content': compute_per_sequence_dinucleotide_content(sequences),
        'Per position nucleotide content': compute_per_position_nucleotide_content(sequences),
        'Per position reversed nucleotide content': compute_per_position_nucleotide_content(sequences, reverse=True),
        'Per sequence GC content': compute_per_sequence_gc_content(sequences),
        'Sequence lenghts': compute_sequence_lengths(sequences),
        'Sequence duplication levels': compute_sequence_duplication_levels(sequences)
    }
    return stats

def evaluate_stats(stats, threshold=0.01):
    """
    Evaluate the given statistics and flag significant differences.
    @param stats: A dictionary containing sequence statistics.
    @param threshold: The p-value threshold for flagging significant differences.
    @return: A dictionary containing evaluated statistics.
                Keys: the same as the input dictionary
                Values: 'passed', 'warning', or 'failed'
    """
    results = {
        'Filename': '',
        'Number of sequences': evaluate_number_of_sequences(stats['Number of sequences']),
        'Number of bases': 'No implemented yet',
        'Per sequence nucleotide content': evaluate_per_sequence_nucleotide_content(stats['Per sequence nucleotide content'], threshold),
        '%GC content': 'No implemented yet',
        'Per sequence dinucleotide content': evaluate_per_sequence_dinucleotide_content(stats['Per sequence dinucleotide content'], threshold),
        'Per position nucleotide content': evaluate_per_position_nucleotide_content(stats['Per position nucleotide content'], threshold),
        'Per position reversed nucleotide content': evaluate_per_position_nucleotide_content(stats['Per position reversed nucleotide content'], threshold),
        'Per sequence GC content': evaluate_per_sequence_gc_content(stats['Per sequence GC content'], threshold),
        'Sequence lenghts': evaluate_sequence_lengths(stats['Sequence lenghts']),
        'Sequence duplication levels': evaluate_sequence_duplication_levels(stats['Sequence duplication levels'])
    }
    return results

def flag_significant_differences(pos_stats, neg_stats, threshold=0.01):
    raise NotImplementedError

def compute_per_sequence_nucleotide_content(sequences):
    """
    Compute the nucleotide content for each sequence in the given list of sequences.
    @return: A dictionary containing the nucleotide content for each sequence.
                Keys: sequence ID
                Values: dictionary {nucleotide: frequency probability (count / len(sequence))}
    """
    nucleotides_per_sequence = {}
    for id, sequence in enumerate(sequences):
        nucleotides_per_sequence[id] = {nucleotide: sequence.count(nucleotide) / len(sequence) for nucleotide in set(sequence)}
    return nucleotides_per_sequence


def compute_per_sequence_dinucleotide_content(sequences):
    """
    Compute the dinucleotide content for each sequence in the given list of sequences.
    @return: A dictionary containing the dinucleotide content for each sequence.
                Keys: sequence ID
                Values: dictionary {dinucleotide: frequency probability (count / len(sequence))}
    """
    dinucleotides_per_sequence = {}
    for id, sequence in enumerate(sequences):
        dinucleotides_per_sequence[id] = {sequence[i:i+2]: sequence.count(sequence[i:i+2]) / len(sequence) for i in range(len(sequence) - 1)}
    return dinucleotides_per_sequence


def compute_per_position_nucleotide_content(sequences, reverse=False):
    """
    Compute the nucleotide content for each position in the given list of sequences. Sequences can be of different lengths.
    @param reverse: If True, compute the nucleotide content for each position in the reversed sequences.
    @return: A dictionary containing the nucleotide content for each position.
                Keys: position
                Values: dictionary {nucleotide: frequency probability (count / len(sequences))}
    """
    nucleotides_per_position = {}
    for sequence in sequences:
        if reverse:
            sequence = sequence[::-1]
        for i, nucleotide in enumerate(sequence):
            if i in nucleotides_per_position:
                nucleotides_per_position[i][nucleotide] = nucleotides_per_position[i].get(nucleotide, 0) + 1
            else:
                nucleotides_per_position[i] = {nucleotide: 1}
    for position in nucleotides_per_position:
        total = sum(nucleotides_per_position[position].values())
        nucleotides_per_position[position] = {nucleotide: count / total for nucleotide, count in nucleotides_per_position[position].items()}
    return nucleotides_per_position


def compute_per_sequence_gc_content(sequences):
    gc_content_per_sequence = {}
    for id, sequence in enumerate(sequences):
        gc_content_per_sequence[id] = (sequence.count('G') + sequence.count('C')) / len(sequence)
    return gc_content_per_sequence

def compute_sequence_lengths(sequences):
    return {id: len(sequence) for id, sequence in enumerate(sequences)}

def compute_sequence_duplication_levels(sequences):
    # TODO: possibly extend to compute the percentage of sequences remaining if deduplicated
    # TODO: from this we can produce a list of overrepresented sequences, plas plot sequence duplication levels vs % total sequences

    sequence_counts = Counter(sequences)
    return dict(sequence_counts)

def evaluate_number_of_sequences(num_sequences):
    if num_sequences < 100:
        return 'failed'
    elif num_sequences < 1000:
        return 'warning'
    else:
        return 'passed'
    
def evaluate_per_sequence_nucleotide_content(nucleotide_content, threshold):
    return 'No implemented yet'

def evaluate_per_sequence_dinucleotide_content(dinucleotide_content, threshold):
    return 'No implemented yet'

def evaluate_per_position_nucleotide_content(nucleotide_content, threshold):
    return 'No implemented yet'

def evaluate_per_sequence_gc_content(gc_content, threshold):
    return 'No implemented yet'

def evaluate_sequence_lengths(sequence_lengths):
    return 'No implemented yet'

def evaluate_sequence_duplication_levels(sequence_duplication_levels):
    return 'No implemented yet'