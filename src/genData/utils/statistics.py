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
        - number of sequences left after deduplication: int
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
        'Number of sequences left after deduplication': len(set(sequences)),
        'Per sequence nucleotide content': compute_per_sequence_nucleotide_content(sequences),
        'Per sequence dinucleotide content': compute_per_sequence_dinucleotide_content(sequences),
        'Per position nucleotide content': compute_per_position_nucleotide_content(sequences),
        'Per position reversed nucleotide content': compute_per_position_nucleotide_content(sequences, reverse=True),
        'Per sequence GC content': compute_per_sequence_gc_content(sequences),
        'Sequence lenghts': compute_sequence_lengths(sequences),
        'Sequence duplication levels': compute_sequence_duplication_levels(sequences)
    }
    stats['Nucleotide content summary'] = get_nucleotide_content_summary(stats['Per sequence nucleotide content'], stats['Unique bases'])
    stats['Dinucleotide content summary'] = get_dinucleotide_content_summary(stats['Per sequence dinucleotide content'], stats['Unique bases'])
    
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

def get_nucleotide_content_summary(content_dict, unique_bases):
    """Helper function to compute summary statistics for nucleotide content"""
    summary = {}
    for nucleotide in unique_bases:
        values = [seq_dict.get(nucleotide, 0) for seq_dict in content_dict.values()]
        values.sort()
        n = len(values)
        summary[nucleotide] = {
            'min': min(values),
            'q1': values[n//4],
            'median': values[n//2],
            'q3': values[n*3//4],
            'max': max(values)
        }
    return summary

def get_dinucleotide_content_summary(content_dict, unique_bases):
    """Helper function to compute summary statistics for dinucleotide content"""
    summary = {}
    dinucleotides = [n1 + n2 for n1 in unique_bases for n2 in unique_bases]
    for dinucleotide in dinucleotides:
        values = [seq_dict.get(dinucleotide, 0) for seq_dict in content_dict.values()]
        values.sort()
        n = len(values)
        summary[dinucleotide] = {
            'min': min(values),
            'q1': values[n//4],
            'median': values[n//2],
            'q3': values[n*3//4],
            'max': max(values)
        }
    return summary



def compute_per_sequence_dinucleotide_content(sequences):
    """
    Compute the dinucleotide content for each sequence in the given list of sequences.
    @return: A dictionary containing the dinucleotide content for each sequence.
                Keys: sequence ID
                Values: dictionary {dinucleotide: frequency probability (count / num_dinucleotides)}
    """
    nucleotides = list(set(''.join(sequences)))
    dinucleotides = [n1 + n2 for n1 in nucleotides for n2 in nucleotides]
    dinucleotides_per_sequence = {}
    for id, sequence in enumerate(sequences):
        dinucleotides_per_sequence[id] = {}
        for i in range(len(sequence) - 1):
            dinucleotide = sequence[i:i + 2]
            dinucleotides_per_sequence[id][dinucleotide] = dinucleotides_per_sequence[id].get(dinucleotide, 0) + 1
        total = sum(dinucleotides_per_sequence[id].values())
        dinucleotides_per_sequence[id] = {dinucleotide: count / total for dinucleotide, count in dinucleotides_per_sequence[id].items()}
        # add zeros for missing dinucleotides
        for dinucleotide in dinucleotides:
            if dinucleotide not in dinucleotides_per_sequence[id]:
                dinucleotides_per_sequence[id][dinucleotide] = 0

    return dinucleotides_per_sequence


def compute_per_position_nucleotide_content(sequences, reverse=False):
    """
    Compute the nucleotide content for each position in the given list of sequences. Sequences can be of different lengths.
    @param reverse: If True, compute the nucleotide content for each position in the reversed sequences.
    @return: A dictionary containing the nucleotide content for each position.
                Keys: position
                Values: dictionary {nucleotide: frequency probability (count / len(sequences))}
    """
    nucleotides = list(set(''.join(sequences)))
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
        # add zeros for missing nucleotides
        for nucleotide in nucleotides:
            if nucleotide not in nucleotides_per_position[position]:
                nucleotides_per_position[position][nucleotide] = 0
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
    # remove sequences that are not duplicated
    sequence_counts = {sequence: count for sequence, count in sequence_counts.items() if count > 1}
    # sort the sequences by their counts
    sequence_counts = dict(sorted(sequence_counts.items(), key=lambda item: item[1], reverse=True))
    
    return sequence_counts