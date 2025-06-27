import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

import pandas as pd
import os
import logging

from genbenchQC.report.sequence_html_report import get_html_template
from genbenchQC.utils.input_utils import read_stats_json, write_stats_json
from genbenchQC.report.plots import (
    plot_plot_basic_descriptive_stats,
    plot_nucleotides, 
    plot_dinucleotides, 
    plot_per_base_sequence_comparison, 
    plot_lengths,
    plot_gc_content,
    plot_duplicates
)

def generate_html_report(stats_dict, output_path):
    """
    Generate an HTML report from the given statistics dictionary.
    Plots are generated using the Plotly library.
    """
    logging.info(f"Generating HTML report: {output_path}")

    # Load the HTML template
    template = get_html_template(stats_dict)

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

    generate_pdf_report(stats1, stats2, results, output_path, threshold=threshold)

def generate_json_report(stats_dict, output_path):
    write_stats_json(stats_dict, output_path)

def generate_text_report(stats_dict, output_path):
    with open(output_path, 'w') as file:
        file.write('Sequence Statistics Report\n')
        for key, value in stats_dict.items():
            file.write(f'\n{key}:\n')
            if isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    file.write(f'  {sub_key}: {sub_value}\n')
            else:
                file.write(f'  {value}\n')

def generate_simple_report(results, output_path):
    with open(output_path, 'w') as file:
        for key, value in results.items():
            if key == 'Filename':
                continue
            _, passed = value
            passed = "Passed" if passed else "Failed"
            file.write(f'{key}: {passed}\n')

def generate_pdf_report(stats1, stats2, results, output_path, threshold):
    logging.info(f"Generating PDF report: {output_path}")

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
        logging.info(f"Duplicate sequences saved to {duplicate_seqs_path}")

    with PdfPages(output_path) as pdf:
        for fig in plots:
            pdf.savefig(fig, bbox_inches='tight')