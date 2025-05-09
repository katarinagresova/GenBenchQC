import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import ranksums, fisher_exact


def flag_significant_differences(sequences1, stats1, sequences2, stats2, threshold=0.01):
    results = {
        'Per sequence nucleotide content': flag_per_sequence_content(stats1, stats2, 'Per sequence nucleotide content', threshold),
        'Per sequence dinucleotide content': flag_per_sequence_content(stats1, stats2, 'Per sequence dinucleotide content', threshold),
        'Per position nucleotide content': flag_per_position_nucleotide_content(stats1, stats2, threshold),
        'Per position reversed nucleotide content': flag_per_position_nucleotide_content(stats1, stats2, threshold, reversed = True),
        'Per sequence GC content': flag_per_sequence_gc_content(stats1, stats2, threshold),
        'Sequence lengths': flag_sequence_lengths(stats1, stats2, threshold),
        'Duplication between positives and negatives': flag_duplication_between_datasets(sequences1, sequences2)
    }

    return results

def flag_per_sequence_content(stats1, stats2, column, threshold, end_position=None):
    
    df1 = pd.DataFrame(stats1[column]).T
    df2 = pd.DataFrame(stats2[column]).T

    if end_position is None:
        end_position = min(len(df1), len(df2))
    
    # get columns names
    bases = list(set(list(df1.columns.values) + list(df2.columns.values)))

    p_values = {}
    for base in bases:

        _, p_value = ranksums(df1[base][:end_position], df2[base][:end_position])
        p_values[base] = p_value

    # Correcting for FDR
    _, new_p_values = fdrcorrection(list(p_values.values()))
    p_values = dict(zip(p_values.keys(), new_p_values))

    # If there is no significant p-value, this test passed
    passed = np.all(np.array(new_p_values) > threshold)
    
    return (p_values, passed)

def flag_per_position_nucleotide_content(stats1, stats2, threshold, reversed=False, end_position=None):
    df1 = pd.DataFrame(stats1['Per position nucleotide content']).T
    df2 = pd.DataFrame(stats2['Per position nucleotide content']).T

    if end_position is None:
        end_position = min(len(df1), len(df2))

    if reversed:
        df1 = pd.DataFrame(stats1['Per position reversed nucleotide content']).T
        df2 = pd.DataFrame(stats2['Per position reversed nucleotide content']).T

    # get columns names
    bases = list(set(list(df1.columns.values) + list(df2.columns.values)))

    p_values = {}
    passed = True
    for base in bases:

        p_values[base] = []
        for i in range(end_position):
            table=[[df1[base][i] * 100, (1 - df1[base][i]) * 100],
                [df2[base][i] * 100, (1 - df2[base][i]) * 100]]

            _, p_value = fisher_exact(table=table) 
            p_values[base].append(p_value)

        # Correcting for FDR per base
        _, p_values[base] = fdrcorrection(p_values[base])

        passed = passed and np.all(np.array(p_values[base]) > threshold)
 
    return (p_values, passed)

    
def flag_per_sequence_gc_content(stats1, stats2, threshold):

    _, p_value = ranksums(list(stats1['Per sequence GC content'].values()), list(stats2['Per sequence GC content'].values()))
    passed = p_value > threshold

    return (p_value, passed)
    
def flag_sequence_lengths(stats1, stats2, threshold):

    _, p_value = ranksums(list(stats1['Sequence lengths'].values()), list(stats2['Sequence lengths'].values()))
    passed = p_value > threshold

    return (p_value, passed)
    
def flag_duplication_between_datasets(sequences1, sequences2):
    passed = len(set(sequences1).intersection(sequences2)) == 0
    return (None, passed)
