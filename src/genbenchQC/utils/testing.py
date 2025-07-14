import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import wasserstein_distance, fisher_exact
import logging


def flag_significant_differences(sequences1, stats1, sequences2, stats2, threshold, end_position=None):
    results = {
        'Unique bases': flag_unique_bases(
            stats1, stats2
        ),
        'Per sequence nucleotide content': flag_per_sequence_content(
            stats1, stats2, 
            column='Per sequence nucleotide content', 
            threshold=threshold
        ),
        'Per sequence dinucleotide content': flag_per_sequence_content(
            stats1, stats2, 
            column='Per sequence dinucleotide content', 
            threshold=threshold)
        ,
        'Per position nucleotide content': flag_per_position_nucleotide_content(
            stats1, stats2, 
            column='Per position nucleotide content', 
            threshold=threshold, 
            end_position=end_position
        ),
        'Per position reversed nucleotide content': flag_per_position_nucleotide_content(
            stats1, stats2, 
            column='Per position reversed nucleotide content', 
            threshold=threshold, 
            end_position=end_position
        ),
        'Per sequence GC content': flag_per_sequence_one_stat(
            stats1, stats2, 
            column='Per sequence GC content', 
            threshold=threshold
        ),
        'Sequence lengths': flag_per_sequence_one_stat(
            stats1, stats2, 
            column='Sequence lengths',
            threshold=threshold
        ),
        'Duplication between labels': flag_duplication_between_datasets(
            sequences1, sequences2
        )
    }

    return results

def flag_unique_bases(stats1, stats2):
    if set(stats1['Unique bases']) == set(stats2['Unique bases']):
        return (None, True)
    else:
        return (None, False)

def flag_per_sequence_content(stats1, stats2, column, threshold):
    
    df1 = stats1[column]
    df2 = stats2[column]
    
    # get columns names
    bases = list(set(list(df1.columns.values) + list(df2.columns.values)))

    distances = {}
    for base in bases:
        if base not in df1 or base not in df2:
            distances[base] = np.inf
        else:
            distances[base] = wasserstein_distance(df1[base], df2[base])
            logging.debug(f"Distance for {base}: {distances[base]} (threshold: {threshold})")
            max_value = max(
                max(df1[base]),
                max(df2[base])
            )
            if max_value > 0:
                distances[base] /= max_value
                logging.debug(f"Max value for {base}: {max_value}")
                logging.debug(f"Distance after normalization for {base}: {distances[base]}")

    passed = np.all(np.array(list(distances.values())) < threshold)
    
    return (distances, passed)

def flag_per_position_nucleotide_content(stats1, stats2, column, threshold, end_position):
    
    df1 = stats1[column]
    df2 = stats2[column]

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

    distance = wasserstein_distance(
        stats1[column].values.flatten(),
        stats2[column].values.flatten()
    )
    logging.debug(f"Distance for {column}: {distance} (threshold: {threshold})")
    max_value = max(
        max(stats1[column].values.flatten()), 
        max(stats2[column].values.flatten())
    )
    if max_value > 0:
        distance /= max_value
        logging.debug(f"Max value for {column}: {max_value}")
        logging.debug(f"Distance after normalization: {distance}")

    passed = distance < threshold

    return (distance, passed)
    
def flag_duplication_between_datasets(sequences1, sequences2):
    duplicates = list(set(sequences1).intersection(sequences2))
    if len(duplicates) > 0:
        return (duplicates, False)
    else:
        return ([], True)
