import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

import pandas as pd

from genbenchQC.report.sequence_html_report import get_html_template

def generate_html_report(stats_dict, output_path):
    """
    Generate an HTML report from the given statistics dictionary.
    Plots are generated using the Plotly library.
    """
    # Load the HTML template
    template = get_html_template(stats_dict)

    with open(output_path, 'w') as file:
        file.write(template)

def generate_dataset_html_report(positive_stats, negative_stats, results, output_path):
    """
    Generate an HTML report comparing the statistics of two datasets.
    Plots are generated using the Plotly library.
    """
    with open(output_path, 'w') as file:
        file.write("Not implemented yet. Producing a PDF report instead.")

    # get the output path without the extension
    output_path = str(output_path).split('.')[0] + '.pdf'

    generate_pdf_report(positive_stats, negative_stats, results, output_path)


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

def generate_pdf_report(pos_stats, neg_stats, results, output_path, threshold=0.01):

    plots = []
    bases = list(set(pos_stats['Unique bases'] + neg_stats['Unique bases']))

    # Plot per sequence nucleotide content
    fig, ax = plt.subplots(nrows=len(bases), ncols=1, figsize=(10, 2*len(bases)))
    plot_composition_comparison(
        pos_stats['Per sequence nucleotide content'],
        neg_stats['Per sequence nucleotide content'],
        bases = bases,
        x_label='Frequency',
        ax=ax,
        stats=results['Per sequence nucleotide content'][0],
        p_value_thresh=threshold
    )
    fig.suptitle('Nucleotide composition')
    plots.append(fig)

    dinucleotides = list(set(list(pos_stats['Per sequence dinucleotide content'][0].keys()) + list(neg_stats['Per sequence dinucleotide content'][0].keys())))
    # Plot per sequence dinucleotide content
    fig, ax = plt.subplots(nrows=len(dinucleotides), ncols=1, figsize=(10, 2*len(dinucleotides)))
    plot_composition_comparison(
        pos_stats['Per sequence dinucleotide content'],
        neg_stats['Per sequence dinucleotide content'],
        bases = dinucleotides,
        x_label='Frequency',
        ax=ax,
        stats=results['Per sequence dinucleotide content'][0],
        p_value_thresh=threshold
    )
    fig.suptitle('Dinucleotide composition')
    plots.append(fig)

    # Plot per position nucleotide content
    fig, ax = plt.subplots(nrows=len(bases), ncols=1, figsize=(10, 2*len(bases)))
    plot_per_base_sequence_comparison(
        pos_stats['Per position nucleotide content'],
        neg_stats['Per position nucleotide content'],
        bases = bases,
        end_position=100,
        x_label='Position in read (bp)',
        ax=ax,
        stats=results['Per position nucleotide content'][0],
        p_value_thresh=threshold
    )
    fig.suptitle('Nucleotide composition per position')
    plots.append(fig)

    # Plot per position reversed nucleotide content
    fig, ax = plt.subplots(nrows=len(bases), ncols=1, figsize=(10, 2*len(bases)))
    plot_per_base_sequence_comparison(
        pos_stats['Per position reversed nucleotide content'],
        neg_stats['Per position reversed nucleotide content'],
        bases = bases,
        end_position=100,
        x_label='Position in read reversed (bp)',
        ax=ax,
        stats=results['Per position reversed nucleotide content'][0],
        p_value_thresh=threshold
    )
    fig.suptitle('Reversed nucleotide composition per position')
    plots.append(fig)

    # Plot length distribution
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 2))
    plot_lenght_comparison(
        pos_stats['Sequence lengths'],
        neg_stats['Sequence lengths'],
        x_label='Length',
        ax=ax,
        stats=results['Sequence lengths'][0],
        p_value_thresh=threshold
    )
    fig.suptitle('Sequence length distribution')
    plots.append(fig)


    with PdfPages(output_path) as pdf:
        for fig in plots:
            pdf.savefig(fig, bbox_inches='tight')

def plot_per_base_sequence_content(df, end_position, title='', ax=None):

    if ax is None:
        fig, ax = plt.subplots(figsize=(15, 5))

    for column in df.columns[:-1]:        
        ax.plot(df.index.to_numpy()[:end_position], df[column].to_numpy()[:end_position] / df['sum'].to_numpy()[:end_position], label=column)
            
    ax.legend()

    ax.set_xlabel('Position in read (bp)')
    ax.set_ylabel('Frequency')
    ax.set_ylim(0, 1)
    ax.set_title(title)

    return ax

def plot_per_base_sequence_comparison(pos_stats, neg_stats, bases, end_position, x_label='', ax=None, stats=None, p_value_thresh=0.01):

    if ax is None:
        fig, ax = plt.subplots(nrows=len(bases), ncols=1, figsize=(15, 3 * len(bases)), sharex=True, sharey=True)

    for index, base in enumerate(bases):

        pos_df = pd.DataFrame(pos_stats).T
        neg_df = pd.DataFrame(neg_stats).T

        ax[index].plot(pos_df.index[:end_position], pos_df[base][:end_position], label=f"positive {base}")
        ax[index].plot(neg_df.index[:end_position], neg_df[base][:end_position], label=f"negative {base}")

        ax[index].set_ylim(-0.1, 1.1)
        ax[index].set_ylabel('Frequency')
        ax[index].legend()
        ax[index].set_title(f'Base: {base}')

        if stats:
            for i in range(end_position):
                
                p_value = stats[base][i]
                if p_value < p_value_thresh:
                    ax[index].plot(i, p_value, 'ro',)

                    # plot salmon background on position with p-value < 0.05
                    ax[index].axvspan(i-0.45, i+0.45, facecolor='salmon', alpha=0.5)

                # add new value to legend
                ax[index].legend([f"positive {base}", f"negative {base}", f"p-value < 0.05"])

    ax[index].set_xlabel(x_label)

    return ax

def plot_composition_comparison(pos_stats, neg_stats, bases, x_label='', ax=None, stats=None, p_value_thresh=0.01):

    pos_df = pd.DataFrame(pos_stats).T
    neg_df = pd.DataFrame(neg_stats).T

    if ax is None:
        fig, ax = plt.subplots(nrows=len(bases), ncols=1, figsize=(10, 2*len(bases)))
    for index, base in enumerate(bases):

        sns.histplot(pos_df[base], ax=ax[index], label='positive', alpha=0.3, stat='probability')
        sns.histplot(neg_df[base], ax=ax[index], label='negative', alpha=0.3, stat='probability')

        if stats:
            p_value = stats.get(base, 0)
            if p_value < p_value_thresh:
                ax[index].text(0.9, 0.01, f"p-value: {p_value:.2e}", ha='center')
                
                # set frame to gray and width to 2
                ax[index].spines['bottom'].set_color('red')
                ax[index].spines['top'].set_color('red')
                ax[index].spines['right'].set_color('red')
                ax[index].spines['left'].set_color('red')
                ax[index].spines['bottom'].set_linewidth(2)
                ax[index].spines['top'].set_linewidth(2)
                ax[index].spines['right'].set_linewidth(2)
                ax[index].spines['left'].set_linewidth(2)

        ax[index].set_xlim(0, 1)
        ax[index].set_ylabel('Seq Frequency')
        ax[index].legend()
        ax[index].set_title(f'Base: {base}')

    ax[index].set_xlabel(x_label)

    return ax


def plot_composition_comparison_boxplot(pos_df, neg_df, title='', x_label='', ax=None, stats=None, p_value_thresh=0.01):

    # get columns names
    bases = pos_df.columns.values
    # remove sum from bases (it doesn't have to be the last value)
    bases = bases[bases != 'sum']

    if ax is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(2 * len(bases), 5))

    new_pos_df = pos_df.assign(label='positive')
    new_neg_df = neg_df.assign(label='negative')

    combined_df = pd.concat([new_pos_df, new_neg_df])
    del new_pos_df, new_neg_df

    for base in bases:
        combined_df[base] = combined_df[base] / combined_df['sum']
        
    combined_df = combined_df.melt(id_vars=['label'], value_vars=bases, var_name='base', value_name='frequency')

    sns.boxplot(data=combined_df, x='base', y='frequency', hue='label', ax=ax)

    for base in bases:

        if stats:
            p_value = stats.get(base, 0)
            if p_value < p_value_thresh:

                # plot p-value on top of the boxplot
                index = list(bases).index(base)
                ax.text(index, -0.15, f"p-value: {p_value:.2e}", ha='center')
                # set background of given column to salmon
                ax.axvspan(index-0.48, index+0.48, facecolor='salmon', alpha=0.5)

    ax.set_title(title)
    ax.set_ylabel('Frequency')
    ax.set_ylim(-0.2, 1.1)
    ax.set_xlabel(x_label)

    return ax

#@DeprecationWarning
def plot_lenght_comparison(pos_stats, neg_stats, title='', x_label='', ax=None, stats=None, p_value_thresh=0.01):

    if ax is None:
        fig, ax = plt.subplots(1, ncols=1, figsize=(10, 2))

    pos_length_counts = pos_stats.values() 
    neg_length_counts = neg_stats.values()

    sns.histplot(pos_length_counts, ax=ax, label='positive', alpha=0.3, stat='probability')
    sns.histplot(neg_length_counts, ax=ax, label='negative', alpha=0.3, stat='probability')

    if stats:
        p_value = stats
        if p_value < p_value_thresh:
            ax.text(0.9, 0.1, f"p-value: {p_value:.2e}", ha='center', transform=ax.transAxes)
            
            # set frame to gray and width to 2
            ax.spines['bottom'].set_color('salmon')
            ax.spines['top'].set_color('salmon')
            ax.spines['right'].set_color('salmon')
            ax.spines['left'].set_color('salmon')
            ax.spines['bottom'].set_linewidth(2)
            ax.spines['top'].set_linewidth(2)
            ax.spines['right'].set_linewidth(2)
            ax.spines['left'].set_linewidth(2)

    #ax.set_xlim(0, 1)
    ax.set_ylabel('Frequency')
    ax.legend()

    ax.set_xlabel(x_label)

    return ax

#@DeprecationWarning
def plot_lenght_comparison_boxplot(pos_df, neg_df, title='', x_label='', ax=None, stats=None, p_value_thresh=0.01):

    if ax is None:
        fig, ax = plt.subplots(1, ncols=1, figsize=(10, 2))

    new_pos_df = pos_df.assign(label='positive')
    new_neg_df = neg_df.assign(label='negative')

    combined_df = pd.concat([new_pos_df, new_neg_df])
    combined_df = combined_df[['sum', 'label']]
    del new_pos_df, new_neg_df

    sns.boxplot(data=combined_df, x='label', y='sum', ax=ax)

    if stats:
        p_value = stats
        if p_value < p_value_thresh:

            ax.text(0.5, 0.1, f"p-value: {p_value:.2e}", ha='center', transform=ax.transAxes)
            
            # set frame to gray and width to 2
            ax.spines['bottom'].set_color('salmon')
            ax.spines['top'].set_color('salmon')
            ax.spines['right'].set_color('salmon')
            ax.spines['left'].set_color('salmon')
            ax.spines['bottom'].set_linewidth(2)
            ax.spines['top'].set_linewidth(2)
            ax.spines['right'].set_linewidth(2)
            ax.spines['left'].set_linewidth(2)

    ax.set_ylabel('Length')

    ax.set_xlabel(x_label)

    return ax
    