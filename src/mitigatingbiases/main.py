from counting import do_counting
from plotting.plotting import plot_per_base_sequence_comparison, plot_composition_comparison, plot_composition_comparison_boxplot, plot_lenght_comparison, plot_lenght_comparison_boxplot
from stats.stats import stats_composition_comparison, stats_per_base_sequence_comparison, stats_lenght_comparison

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Bio import SeqIO
import argparse

def parse_fasta(filename):
    """
    Parses a fasta file and returns a generator of SeqIO objects
    """

    with open(filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            yield record.seq.upper()

def print_passed(stats):
    # Print stats. Failed in red and passed in green
    for key, value in stats.items():
        if value:
            print(f'{key}: \033[92mPASSED\033[0m')
        else:
            print(f'{key}: \033[91mFAILED\033[0m')

def save_plots(plots, output_file):
    with PdfPages(output_file) as pdf:
        for fig in plots:
            pdf.savefig(fig, bbox_inches='tight')

def parse_args():
    parser = argparse.ArgumentParser(description='Mitigating biases')
    parser.add_argument('--pos', type=str, required=True, help='Path to positive sequences')
    parser.add_argument('--neg', type=str, required=True, help='Path to negative sequences')
    parser.add_argument('--output', type=str, required=True, help='Path to output PDF report file')
    parser.add_argument('--p_value_thresh', type=float, default=0.01, help='P-value threshold')
    return parser.parse_args()

def get_plots(
    pos_nucleotides_df,
    neg_nucleotides_df,
    pos_dinucleotides_df,
    neg_dinucleotides_df,
    pos_nucleotides_per_position_df,
    neg_nucleotides_per_position_df,
    pos_nucleotides_per_position_reversed_df,
    neg_nucleotides_per_position_reversed_df,
    p_value_thresh=0.01,
    stats=None
):
    plots = []

    # get bcolumns names
    bases = pos_nucleotides_df.columns.values
    # remove sum from bases (it doesn't have to be the last value)
    bases = bases[bases != 'sum']

    # Plot nucleotides
    fig, ax = plt.subplots(nrows=len(bases), ncols=1, figsize=(10, 2*len(bases)))
    plot_composition_comparison(
        pos_nucleotides_df,
        neg_nucleotides_df,
        x_label='Frequency',
        ax=ax,
        stats=stats['Nucleotide composition'],
        p_value_thresh=p_value_thresh)
    fig.suptitle('Nucleotide composition')
    plots.append(fig)

    # get bcolumns names
    bases = pos_dinucleotides_df.columns.values
    # remove sum from bases (it doesn't have to be the last value)
    bases = bases[bases != 'sum']

    # Plot dinucleotides
    fig, ax = plt.subplots(nrows=len(bases), ncols=1, figsize=(10, 2*len(bases)))
    plot_composition_comparison(
        pos_dinucleotides_df,
        neg_dinucleotides_df,
        x_label='Frequency',
        ax=ax,
        stats=stats['Dinucleotide composition'],
        p_value_thresh=p_value_thresh)
    fig.suptitle('Dinucleotide composition')
    plots.append(fig)

    # get columns names
    bases = pos_nucleotides_per_position_df.columns.values
    # remove sum from bases (it doesn't have to be the last value)
    bases = bases[bases != 'sum']

    # Plot nucleotides per position
    fig, ax = plt.subplots(nrows=len(bases), ncols=1, figsize=(10, 2*len(bases)))
    plot_per_base_sequence_comparison(
        pos_nucleotides_per_position_df, 
        neg_nucleotides_per_position_df, 
        end_position=100,  
        x_label='Position in read (bp)',
        ax=ax,
        stats=stats['Nucleotide per position'],
        p_value_thresh=p_value_thresh)
    fig.suptitle('Nucleotide per position')
    plots.append(fig)

    # get columns names
    bases = pos_nucleotides_per_position_reversed_df.columns.values
    # remove sum from bases (it doesn't have to be the last value)
    bases = bases[bases != 'sum']

    # Plot nucleotides per position reversed
    fig, ax = plt.subplots(nrows=len(bases), ncols=1, figsize=(10, 2 * len(bases)))
    plot_per_base_sequence_comparison(
        pos_nucleotides_per_position_reversed_df, 
        neg_nucleotides_per_position_reversed_df, 
        end_position=100, 
        x_label='Position in read (bp)',
        ax=ax,
        stats=stats['Nucleotide per position reversed'],
        p_value_thresh=p_value_thresh)
    fig.suptitle('Nucleotide per position reversed')
    plots.append(fig)

    # Plot length distribution
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 3))
    plot_lenght_comparison(
        pos_nucleotides_df, 
        neg_nucleotides_df, 
        x_label='Sequence length', 
        ax=ax,
        stats=stats['Length distribution'],
        p_value_thresh=p_value_thresh)
    fig.suptitle('Length distribution')
    plots.append(fig)

    return plots

def get_box_plots(
    pos_nucleotides_df,
    neg_nucleotides_df,
    pos_dinucleotides_df,
    neg_dinucleotides_df,
    pos_nucleotides_per_position_df,
    neg_nucleotides_per_position_df,
    pos_nucleotides_per_position_reversed_df,
    neg_nucleotides_per_position_reversed_df,
    p_value_thresh=0.01,
    stats=None
):
    plots = []
    
    # get bcolumns names
    bases = pos_nucleotides_df.columns.values
    # remove sum from bases (it doesn't have to be the last value)
    bases = bases[bases != 'sum']

    # Plot nucleotides
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(2*len(bases), 6))
    plot_composition_comparison_boxplot(
        pos_nucleotides_df,
        neg_nucleotides_df,
        x_label='Nucleotides',
        ax=ax,
        stats=stats['Nucleotide composition'],
        p_value_thresh=p_value_thresh
    )
    fig.suptitle('Nucleotide composition')
    plots.append(fig)

    # get bcolumns names
    bases = pos_dinucleotides_df.columns.values
    # remove sum from bases (it doesn't have to be the last value)
    bases = bases[bases != 'sum']

    # Plot dinucleotides
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(2*len(bases), 6))
    plot_composition_comparison_boxplot(
        pos_dinucleotides_df,
        neg_dinucleotides_df,
        x_label='Dinucleotides',
        ax=ax,
        stats=stats['Dinucleotide composition'],
        p_value_thresh=p_value_thresh
    )
    fig.suptitle('Dinucleotide composition')
    plots.append(fig)

    # get columns names
    bases = pos_nucleotides_per_position_df.columns.values
    # remove sum from bases (it doesn't have to be the last value)
    bases = bases[bases != 'sum']

    # Plot nucleotides per position
    fig, ax = plt.subplots(nrows=len(bases), ncols=1, figsize=(10, 2*len(bases)))
    plot_per_base_sequence_comparison(
        pos_nucleotides_per_position_df, 
        neg_nucleotides_per_position_df, 
        end_position=100,  
        x_label='Position in read (bp)',
        ax=ax,
        stats=stats['Nucleotide per position'],
        p_value_thresh=p_value_thresh)
    fig.suptitle('Nucleotide per position')
    plots.append(fig)

    # get columns names
    bases = pos_nucleotides_per_position_reversed_df.columns.values
    # remove sum from bases (it doesn't have to be the last value)
    bases = bases[bases != 'sum']

    # Plot nucleotides per position reversed
    fig, ax = plt.subplots(nrows=len(bases), ncols=1, figsize=(10, 2 * len(bases)))
    plot_per_base_sequence_comparison(
        pos_nucleotides_per_position_reversed_df, 
        neg_nucleotides_per_position_reversed_df, 
        end_position=100, 
        x_label='Position in read (bp)',
        ax=ax,
        stats=stats['Nucleotide per position reversed'],
        p_value_thresh=p_value_thresh)
    fig.suptitle('Nucleotide per position reversed')
    plots.append(fig)

    # Plot length distribution
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    plot_lenght_comparison_boxplot(
        pos_nucleotides_df, 
        neg_nucleotides_df, 
        x_label='Dataset', 
        ax=ax,
        stats=stats['Length distribution'],
        p_value_thresh=p_value_thresh)
    fig.suptitle('Length distribution')
    plots.append(fig)

    return plots


def get_stats(
    pos_nucleotides_df,
    neg_nucleotides_df,
    pos_dinucleotides_df,
    neg_dinucleotides_df,
    pos_nucleotides_per_position_df,
    neg_nucleotides_per_position_df,
    pos_nucleotides_per_position_reversed_df,
    neg_nucleotides_per_position_reversed_df,
    p_value_thresh=0.01
):
    stats = {}
    passed = {}

    stats['Nucleotide composition'], passed['Nucleotide composition'] = stats_composition_comparison(
        pos_nucleotides_df,
        neg_nucleotides_df,
        end_position=100,
        p_value_thresh=p_value_thresh
    )
    stats['Dinucleotide composition'], passed['Dinucleotide composition'] = stats_composition_comparison(
        pos_dinucleotides_df,
        neg_dinucleotides_df,
        end_position=100,
        p_value_thresh=p_value_thresh
    )
    stats['Nucleotide per position'], passed['Nucleotide per position'] = stats_per_base_sequence_comparison(
        pos_nucleotides_per_position_df,
        neg_nucleotides_per_position_df,
        end_position=100,
        p_value_thresh=p_value_thresh
    )
    stats['Nucleotide per position reversed'], passed['Nucleotide per position reversed'] = stats_per_base_sequence_comparison(
        pos_nucleotides_per_position_reversed_df,
        neg_nucleotides_per_position_reversed_df,
        end_position=100,
        p_value_thresh=p_value_thresh
    )
    stats['Length distribution'], passed['Length distribution'] = stats_lenght_comparison(
        pos_nucleotides_df['sum'],
        neg_nucleotides_df['sum'],
        p_value_thresh=p_value_thresh
    )

    return stats, passed

def main():

    args = parse_args()

    pos_seq = parse_fasta(args.pos)
    neg_seq = parse_fasta(args.neg)

    pos_nucleotides_df, pos_dinucleotides_df, pos_nucleotides_per_position_df, pos_nucleotides_per_position_reversed_df = do_counting(pos_seq)
    neg_nucleotides_df, neg_dinucleotides_df, neg_nucleotides_per_position_df, neg_nucleotides_per_position_reversed_df = do_counting(neg_seq)

    stats, passed = get_stats(
        pos_nucleotides_df,
        neg_nucleotides_df,
        pos_dinucleotides_df,
        neg_dinucleotides_df,
        pos_nucleotides_per_position_df,
        neg_nucleotides_per_position_df,
        pos_nucleotides_per_position_reversed_df,
        neg_nucleotides_per_position_reversed_df,
        p_value_thresh=args.p_value_thresh
    )

    plots = get_box_plots(
        pos_nucleotides_df,
        neg_nucleotides_df,
        pos_dinucleotides_df,
        neg_dinucleotides_df,
        pos_nucleotides_per_position_df,
        neg_nucleotides_per_position_df,
        pos_nucleotides_per_position_reversed_df,
        neg_nucleotides_per_position_reversed_df,
        stats=stats,
        p_value_thresh=args.p_value_thresh
    )

    print_passed(passed)
    save_plots(plots, output_file=args.output)

if __name__ == "__main__":
    main()