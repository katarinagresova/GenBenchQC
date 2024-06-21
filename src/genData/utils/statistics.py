import numpy as np
import pandas as pd
from collections import Counter

from genData.utils.fasta_utils import read_fasta

def compute_sequence_statistics(fasta_file):
    """
    Compute various statistics from the given list of sequences.
    @param sequences: A list of sequences.
    @return: A dictionary containing the following statistics:
        - Filename: str
        - Number of sequences: int
        - Number of bases: int
        - Unique bases: list of str
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
        'Unique bases': list(set(''.join(sequences))),
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

def compute_per_sequence_nucleotide_content(sequences):
    """
    Compute the nucleotide content for each sequence in the given list of sequences.
    @return: A dictionary containing the nucleotide content for each sequence.
                Keys: sequence ID
                Values: dictionary {nucleotide: frequency probability (count / len(sequence))}
    """
    nucleotides = list(set(''.join(sequences)))
    nucleotides_per_sequence = {}
    for id, sequence in enumerate(sequences):
        nucleotides_per_sequence[id] = {nucleotide: sequence.count(nucleotide) / len(sequence) for nucleotide in set(sequence)}
        # add zeros for missing nucleotides
        for nucleotide in nucleotides:
            if nucleotide not in nucleotides_per_sequence[id]:
                nucleotides_per_sequence[id][nucleotide] = 0
    return nucleotides_per_sequence


def compute_per_sequence_dinucleotide_content(sequences):
    """
    Compute the dinucleotide content for each sequence in the given list of sequences.
    @return: A dictionary containing the dinucleotide content for each sequence.
                Keys: sequence ID
                Values: dictionary {dinucleotide: frequency probability (count / num_dinucleotides)}
    """
    dinucleotides_per_sequence = {}
    for id, sequence in enumerate(sequences):
        dinucleotides_per_sequence[id] = {}
        for i in range(len(sequence) - 1):
            dinucleotide = sequence[i:i + 2]
            dinucleotides_per_sequence[id][dinucleotide] = dinucleotides_per_sequence[id].get(dinucleotide, 0) + 1
        total = sum(dinucleotides_per_sequence[id].values())
        dinucleotides_per_sequence[id] = {dinucleotide: count / total for dinucleotide, count in dinucleotides_per_sequence[id].items()}
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
    """
    Compute the duplication levels for each sequence in the given list of sequences.
    @param sequences: A list of sequences.
    @return: A dictionary containing the duplication levels for duplicated sequences. Unique sequences are not included.
    """
    # TODO: possibly extend to compute the percentage of sequences remaining if deduplicated
    # TODO: from this we can produce a list of overrepresented sequences, plus plot sequence duplication levels vs % total sequences

    sequence_counts = Counter(sequences)
    return {sequence: count for sequence, count in sequence_counts.items() if count > 1}