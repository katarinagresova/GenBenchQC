import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

import pandas as pd
import os
from pathlib import Path

from genbenchQC.report.sequence_html_report import get_html_template
from genbenchQC.utils.input_utils import write_stats_json
from genbenchQC.report.dataset_plots import (
    plot_plot_basic_descriptive_stats,
    plot_nucleotides, 
    plot_dinucleotides, 
    plot_per_base_sequence_comparison, 
    plot_lengths,
    plot_gc_content,
    plot_duplicates
)

from genbenchQC.report.sequences_plots import (
    hist_plot_one_stat,
)

def generate_plots(stats_dict, output_path):
    """
    Generate a plots from the given statistics dictionary.
    """

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

    return plots_paths

def generate_html_report(stats_dict, output_path, plots_path):
    """
    Generate an HTML report from the given statistics dictionary.
    Plots are generated using the Plotly library.
    """

    plots_path.mkdir(parents=True, exist_ok=True)

    # generate 
    plots_paths = generate_plots(stats_dict, plots_path)

    # Load the HTML template
    template = get_html_template(stats_dict, plots_paths)

    with open(output_path, 'w') as file:
        file.write(template)

def generate_dataset_html_report(stats1, stats2, results, output_path, threshold):
    """
    Generate an HTML report comparing the statistics of two datasets.
    Plots are generated using the Plotly library.
    """
    with open(output_path, 'w') as file:
        file.write("Not implemented yet. Producing a PDF report instead.")

    # get the output path without the extension
    output_path = os.path.splitext(output_path)[0] + '.pdf'

    generate_dataset_pdf_report(stats1, stats2, results, output_path, threshold=threshold)

def generate_json_report(stats_dict, output_path):
    write_stats_json(stats_dict, output_path)

def generate_simple_report(results, output_path):

    # construct table from results - Name and Passed/Failed status
    df = pd.DataFrame({
        'Name': [key for key in results.keys()],
        'Passed': [value[1] for value in results.values()]
    })
    # convert boolean to string
    df['Passed'] = df['Passed'].apply(lambda x: 'Passed' if x else 'Failed')
    # save to csv
    df.to_csv(output_path, index=False, header=False)

def generate_dataset_pdf_report(stats1, stats2, results, output_path, threshold):

    plots = []
    bases_overlap = list(set(stats1.stats['Unique bases']) & set(stats2.stats['Unique bases']))

    # Plot basic descriptive statistics
    fig = plot_plot_basic_descriptive_stats(stats1, stats2)
    plots.append(fig)
    plt.close(fig)

    # Plot per sequence nucleotide content
    fig = plot_nucleotides(
        stats1,
        stats2,
        nucleotides = bases_overlap,
        result=results['Per sequence nucleotide content'],
        dist_thresh=threshold,
        plot_type='violin'
    )
    plots.append(fig)
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
    plots.append(fig)
    plt.close(fig)
    
    # Plot per position nucleotide content
    fig = plot_per_base_sequence_comparison(
        stats1,
        stats2,
        stats_name='Per position nucleotide content',
        result=results['Per position nucleotide content'],
        p_value_thresh=threshold,
        nucleotides = bases_overlap,
        end_position=100,
        x_label='Position in read (bp)',
        title='Nucleotide composition per position',
    )
    plots.append(fig)
    plt.close(fig)

    # Plot per reversed position nucleotide content
    fig = plot_per_base_sequence_comparison(
        stats1,
        stats2,
        stats_name='Per position reversed nucleotide content',
        result=results['Per position reversed nucleotide content'],
        p_value_thresh=threshold,
        nucleotides = bases_overlap,
        end_position=100,
        x_label='Position in read reversed (bp)',
        title='Reversed nucleotide composition per position',
    )
    plots.append(fig)
    plt.close(fig)

    # Plot length distribution
    fig = plot_lengths(
        stats1,
        stats2,
        result=results['Sequence lengths'],
        dist_thresh=threshold,
    )
    plots.append(fig)
    plt.close(fig)

    # Plot per sequence GC content
    fig = plot_gc_content(
        stats1,
        stats2,
        result=results['Per sequence GC content'],
        dist_thresh=threshold,
    )
    plots.append(fig)
    plt.close(fig)

    fig = plot_duplicates(result=results['Duplication between labels'])
    if fig:
        plots.append(fig)
        plt.close(fig)

        # save duplicate sequences to a file
        duplicate_seqs = results['Duplication between labels'][0]
        # remove extension from output path, add '_duplicates.txt'
        duplicate_seqs_path = os.path.splitext(output_path)[0] + '_duplicates.txt'
        with open(duplicate_seqs_path, 'w') as f:
            for seq in duplicate_seqs:
                f.write(f"{seq}\n")
        print(f"Duplicate sequences saved to {duplicate_seqs_path}")

    with PdfPages(output_path) as pdf:
        for fig in plots:
            pdf.savefig(fig, bbox_inches='tight')