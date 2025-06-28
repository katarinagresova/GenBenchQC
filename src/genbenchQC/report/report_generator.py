import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

import pandas as pd
import os
from pathlib import Path
import logging

from genbenchQC.report.sequence_html_report import get_sequence_html_template
from genbenchQC.report.dataset_html_report import get_dataset_html_template
from genbenchQC.utils.input_utils import write_stats_json
from genbenchQC.report.dataset_plots import (
    plot_nucleotides, 
    plot_dinucleotides, 
    plot_per_base_sequence_comparison, 
    plot_lengths,
    plot_gc_content,
)

from genbenchQC.report.sequences_plots import (
    hist_plot_one_stat,
    violin_plot_nucleotides,
    violin_plot_dinucleotides,
    plot_per_position_nucleotide_content
)

def generate_sequence_plots(stats_dict, output_path, end_position):
    """
    Generate a plots from the given statistics dictionary.
    """

    logging.info(f"Generating PNG plots at: {output_path}")

    nucleotides = sorted(stats_dict['Unique bases'])

    plots_paths = {}
    fig = hist_plot_one_stat(
        stats_dict,
        stats_name='Sequence lengths',
        x_label='Sequence lengths',
        title='Sequence lengths'
    )

    plots_paths['Sequence lengths'] = Path(output_path.name) / 'sequence_lengths.png'
    fig.savefig(output_path / 'sequence_lengths.png', bbox_inches='tight')
    plt.close(fig)

    fig = hist_plot_one_stat(
        stats_dict,
        stats_name='Per sequence GC content',
        x_label='GC content (%)',
        title='Per sequence GC content'
    )
    plots_paths['Per sequence GC content'] = Path(output_path.name) / 'per_sequence_gc_content.png'
    fig.savefig(output_path / 'per_sequence_gc_content.png', bbox_inches='tight')
    plt.close(fig)

    fig = violin_plot_nucleotides(stats_dict)
    plots_paths['Per sequence nucleotide content'] = Path(output_path.name) / 'per_sequence_nucleotide_content.png'
    fig.savefig(output_path / 'per_sequence_nucleotide_content.png', bbox_inches='tight')
    plt.close(fig)  

    fig = violin_plot_dinucleotides(stats_dict, nucleotides=nucleotides)
    plots_paths['Per sequence dinucleotide content'] = Path(output_path.name) / 'per_sequence_dinucleotide_content.png'
    fig.savefig(output_path / 'per_sequence_dinucleotide_content.png', bbox_inches='tight')
    plt.close(fig)

    fig = plot_per_position_nucleotide_content(
        stats_dict,
        stat_name='Per position nucleotide content',
        nucleotides=nucleotides,
        end_position=end_position,
    )
    plots_paths['Per position nucleotide content'] = Path(output_path.name) / 'per_position_nucleotide_content.png'
    fig.savefig(output_path / 'per_position_nucleotide_content.png', bbox_inches='tight')
    plt.close(fig)

    fig = plot_per_position_nucleotide_content(
        stats_dict,
        stat_name='Per position reversed nucleotide content',
        nucleotides=nucleotides,
        end_position=end_position,
    )
    plots_paths['Per position reversed nucleotide content'] = Path(output_path.name) / 'per_position_reversed_nucleotide_content.png'
    fig.savefig(output_path / 'per_position_reversed_nucleotide_content.png', bbox_inches='tight')
    plt.close(fig)

    return plots_paths

def generate_sequence_html_report(stats_dict, output_path, plots_path, end_position):
    """
    Generate an HTML report from the given statistics dictionary.
    Plots are generated using the Plotly library.
    """
    logging.info(f"Generating HTML report: {output_path}")

    plots_path.mkdir(parents=True, exist_ok=True)

    # generate 
    plots_paths = generate_sequence_plots(stats_dict, plots_path, end_position=end_position)

    # Load the HTML template
    template = get_sequence_html_template(stats_dict, plots_paths)

    with open(output_path, 'w') as file:
        file.write(template)

def generate_dataset_html_report(stats1, stats2, results, output_path, plots_path, threshold, end_position):
    """
    Generate an HTML report comparing the statistics of two datasets.
    """
    plots_path.mkdir(parents=True, exist_ok=True)

    # generate 
    plots_paths = generate_dataset_plots(stats1, stats2, results, plots_path, threshold, end_position)

    # Load the HTML template
    template = get_dataset_html_template(stats1, stats2, plots_paths, results)

    with open(output_path, 'w') as file:
        file.write(template)

    if results['Duplication between labels'][0]:
        # save duplicate sequences to a file
        duplicate_seqs = results['Duplication between labels'][0]
        # remove extension from output path, add '_duplicates.txt'
        duplicate_seqs_path = os.path.splitext(output_path)[0] + '_duplicates.txt'
        with open(duplicate_seqs_path, 'w') as f:
            for seq in duplicate_seqs:
                f.write(f"{seq}\n")
        logging.info(f"Duplicate sequences saved to {duplicate_seqs_path}")

def generate_json_report(stats_dict, output_path):
    write_stats_json(stats_dict, output_path)

def generate_simple_report(results, output_path):

    logging.info(f"Generating simple report: {output_path}")

    # construct table from results - Name and Passed/Failed status
    df = pd.DataFrame({
        'Name': [key for key in results.keys()],
        'Passed': [value[1] for value in results.values()]
    })
    # convert boolean to string
    df['Passed'] = df['Passed'].apply(lambda x: 'Passed' if x else 'Failed')
    # save to csv
    df.to_csv(output_path, index=False, header=False)

def generate_dataset_plots(stats1, stats2, results, output_path, threshold, end_position):

    logging.info(f"Generating PNG plots at: {output_path}")

    plots_paths = {}

    bases_overlap = list(set(stats1.stats['Unique bases']) & set(stats2.stats['Unique bases']))

    # Plot per sequence nucleotide content
    fig = plot_nucleotides(
        stats1,
        stats2,
        nucleotides = bases_overlap,
        result=results['Per sequence nucleotide content'],
        dist_thresh=threshold,
        plot_type='violin'
    )
    plots_paths['Per sequence nucleotide content'] = Path(output_path.name) / 'per_sequence_nucleotide_content.png'
    fig.savefig(output_path / 'per_sequence_nucleotide_content.png', bbox_inches='tight')
    plt.close(fig)

    # Plot per sequence dinucleotide content
    fig = plot_dinucleotides(
        stats1,
        stats2,
        nucleotides = bases_overlap,
        result=results['Per sequence dinucleotide content'],
        dist_thresh=threshold,
        plot_type='violin'
    )
    plots_paths['Per sequence dinucleotide content'] = Path(output_path.name) / 'per_sequence_dinucleotide_content.png'
    fig.savefig(output_path / 'per_sequence_dinucleotide_content.png', bbox_inches='tight')
    plt.close(fig)
    
    # Plot per position nucleotide content
    fig = plot_per_base_sequence_comparison(
        stats1,
        stats2,
        stats_name='Per position nucleotide content',
        result=results['Per position nucleotide content'],
        p_value_thresh=threshold,
        nucleotides = bases_overlap,
        end_position=end_position,
        x_label='Position in sequence',
        title='Nucleotide composition per position',
    )
    plots_paths['Per position nucleotide content'] = Path(output_path.name) / 'per_position_nucleotide_content.png'
    fig.savefig(output_path / 'per_position_nucleotide_content.png', bbox_inches='tight')   
    plt.close(fig)

    # Plot per reversed position nucleotide content
    fig = plot_per_base_sequence_comparison(
        stats1,
        stats2,
        stats_name='Per position reversed nucleotide content',
        result=results['Per position reversed nucleotide content'],
        p_value_thresh=threshold,
        nucleotides = bases_overlap,
        end_position=end_position,
        x_label='Position in reversed sequence',
        title='Reversed nucleotide composition per position',
    )
    plots_paths['Per position reversed nucleotide content'] = Path(output_path.name) / 'per_position_reversed_nucleotide_content.png'
    fig.savefig(output_path / 'per_position_reversed_nucleotide_content.png', bbox_inches='tight')
    plt.close(fig)

    # Plot length distribution
    fig = plot_lengths(
        stats1,
        stats2,
        result=results['Sequence lengths'],
        dist_thresh=threshold,
    )
    plots_paths['Sequence lengths'] = Path(output_path.name) / 'sequence_lengths.png'
    fig.savefig(output_path / 'sequence_lengths.png', bbox_inches='tight')
    plt.close(fig)

    # Plot per sequence GC content
    fig = plot_gc_content(
        stats1,
        stats2,
        result=results['Per sequence GC content'],
        dist_thresh=threshold,
    )
    plots_paths['Per sequence GC content'] = Path(output_path.name) / 'per_sequence_gc_content.png'
    fig.savefig(output_path / 'per_sequence_gc_content.png', bbox_inches='tight')
    plt.close(fig)

    return plots_paths