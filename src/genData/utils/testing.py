import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import ranksums, fisher_exact

from genData.utils.fasta_utils import read_fasta


def flag_significant_differences(pos_fasta, pos_stats, neg_fasta, neg_stats, threshold=0.01):
    results = {
        'Per sequence nucleotide content': flag_per_sequence_content(pos_stats, neg_stats, 'Per sequence nucleotide content', threshold),
        'Per sequence dinucleotide content': flag_per_sequence_content(pos_stats, neg_stats, 'Per sequence dinucleotide content', threshold),
        'Per position nucleotide content': flag_per_position_nucleotide_content(pos_stats, neg_stats, threshold),
        'Per position reversed nucleotide content': flag_per_position_nucleotide_content(pos_stats, neg_stats, threshold, reversed = True),
        'Per sequence GC content': flag_per_sequence_gc_content(pos_stats, neg_stats, threshold),
        'Sequence lengths': flag_sequence_lengths(pos_stats, neg_stats, threshold),
        'Duplication between positives and negatives': flag_duplication_between_datasets(pos_fasta, neg_fasta)
    }

    return results

def flag_per_sequence_content(pos_stats, neg_stats, column, threshold, end_position=None):
    
    pos_df = pd.DataFrame(pos_stats[column]).T
    neg_df = pd.DataFrame(neg_stats[column]).T

    if end_position is None:
        end_position = min(len(pos_df), len(neg_df))
    
    # get columns names
    bases = list(set(list(pos_df.columns.values) + list(neg_df.columns.values)))

    p_values = {}
    for base in bases:

        _, p_value = ranksums(pos_df[base][:end_position], neg_df[base][:end_position])
        p_values[base] = p_value

    # Correcting for FDR
    _, new_p_values = fdrcorrection(list(p_values.values()))
    p_values = dict(zip(p_values.keys(), new_p_values))

    # If there is no significant p-value, this test passed
    passed = np.all(np.array(new_p_values) > threshold)
    
    return (p_values, passed)

def flag_per_position_nucleotide_content(pos_stats, neg_stats, threshold, reversed=False, end_position=None):
    pos_df = pd.DataFrame(pos_stats['Per position nucleotide content']).T
    neg_df = pd.DataFrame(neg_stats['Per position nucleotide content']).T

    if end_position is None:
        end_position = min(len(pos_df), len(neg_df))

    if reversed:
        pos_df = pd.DataFrame(pos_stats['Per position reversed nucleotide content']).T
        neg_df = pd.DataFrame(neg_stats['Per position reversed nucleotide content']).T

    # get columns names
    bases = list(set(list(pos_df.columns.values) + list(neg_df.columns.values)))

    p_values = {}
    passed = True
    for base in bases:

        p_values[base] = []
        for i in range(end_position):
            table=[[pos_df[base][i] * 100, (1 - pos_df[base][i]) * 100],
                [neg_df[base][i] * 100, (1 - neg_df[base][i]) * 100]]

            _, p_value = fisher_exact(table=table) 
            p_values[base].append(p_value)

        # Correcting for FDR per base
        _, p_values[base] = fdrcorrection(p_values[base])

        passed = passed and np.all(np.array(p_values[base]) > threshold)
 
    return (p_values, passed)

    
def flag_per_sequence_gc_content(pos_stats, neg_stats, threshold):

    _, p_value = ranksums(list(pos_stats['Per sequence GC content'].values()), list(neg_stats['Per sequence GC content'].values()))
    passed = p_value > threshold

    return (p_value, passed)
    
def flag_sequence_lengths(pos_stats, neg_stats, threshold):

    _, p_value = ranksums(list(pos_stats['Sequence lengths'].values()), list(neg_stats['Sequence lengths'].values()))
    passed = p_value > threshold

    return (p_value, passed)
    
def flag_duplication_between_datasets(pos_fasta, neg_fasta):
    pos_sequences = read_fasta(pos_fasta)
    neg_sequences = read_fasta(neg_fasta)

    passed = len(set(pos_sequences).intersection(neg_sequences)) == 0

    return (None, passed)
