import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import wasserstein_distance, fisher_exact


def flag_significant_differences(sequences1, stats1, sequences2, stats2, threshold):
    results = {
        'Unique bases': flag_unique_bases(stats1, stats2),
        'Per sequence nucleotide content': flag_per_sequence_content(stats1, stats2, 'Per sequence nucleotide content', threshold=threshold),
        'Per sequence dinucleotide content': flag_per_sequence_content(stats1, stats2, 'Per sequence dinucleotide content', threshold=threshold),
        'Per position nucleotide content': flag_per_position_nucleotide_content(stats1, stats2, column='Per position nucleotide content', threshold=threshold),
        'Per position reversed nucleotide content': flag_per_position_nucleotide_content(stats1, stats2, column='Per position reversed nucleotide content', threshold=threshold),
        'Per sequence GC content': flag_per_sequence_one_stat(stats1, stats2, column='Per sequence GC content', threshold=threshold),
        'Sequence lengths': flag_per_sequence_one_stat(stats1, stats2, column='Sequence lengths', threshold=threshold),
        'Duplication between labels': flag_duplication_between_datasets(sequences1, sequences2)
    }

    return results

def flag_unique_bases(stats1, stats2):
    if set(stats1['Unique bases']) == set(stats2['Unique bases']):
        return (None, True)
    else:
        return (None, False)

def flag_per_sequence_content(stats1, stats2, column, threshold, end_position=None):
    
    df1 = pd.DataFrame(stats1[column]).T
    df2 = pd.DataFrame(stats2[column]).T

    if end_position is None:
        end_position = min(len(df1), len(df2))
    
    # get columns names
    bases = list(set(list(df1.columns.values) + list(df2.columns.values)))

    distances = {}
    for base in bases:
        if base not in df1 or base not in df2:
            distances[base] = np.inf
        else:
            df1_base = df1[base][:end_position]
            df2_base = df2[base][:end_position]
            distances[base] = wasserstein_distance(df1_base, df2_base)

    passed = np.all(np.array(list(distances.values())) < threshold)
    
    return (distances, passed)

def flag_per_position_nucleotide_content(stats1, stats2, column, threshold, end_position=None):
    
    df1 = pd.DataFrame(stats1[column]).T
    df2 = pd.DataFrame(stats2[column]).T

    if end_position is None:
        end_position = min(len(df1), len(df2))

    # get second end position - where one of the stats contains less then 50% values
    lengths = np.array(list(stats1['Sequence lengths'].values()))
    lengths_75th = np.percentile(lengths, 75)
    # round to nearest integer
    lengths_75th = int(np.round(lengths_75th))
    
    end_position = min(end_position, lengths_75th)

    # get columns names
    bases = list(set(list(df1.columns.values) + list(df2.columns.values)))

    p_values = {}
    passed = True
    for base in bases:

        p_values[base] = []
        for i in range(end_position):
            if base not in df1 or base not in df2:
                p_values[base].append(np.inf)
            else:
                df1_base = df1[base][i]
                df2_base = df2[base][i]
                table=[[df1_base * 100, (1 - df1_base) * 100],
                    [df2_base * 100, (1 - df2_base) * 100]]

                _, p_value = fisher_exact(table=table) 
                p_values[base].append(p_value)

        # Correcting for FDR per base
        _, p_values[base] = fdrcorrection(p_values[base])

        passed = passed and np.all(np.array(p_values[base]) > threshold)
 
    return (p_values, passed)

    
def flag_per_sequence_one_stat(stats1, stats2, column, threshold):

    distance = wasserstein_distance(list(stats1[column].values()), list(stats2[column].values()))
    passed = distance < threshold

    return (distance, passed)
    
def flag_duplication_between_datasets(sequences1, sequences2):
    diplicates = set(sequences1).intersection(sequences2)
    if len(diplicates) > 0:
        return (diplicates, False)
    else:
        return (None, True)
