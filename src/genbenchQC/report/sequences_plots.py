import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging

def hist_plot_one_stat(stats, stats_name, x_label='', title=''):
    """
    Plot a single statistic from two stats objects.
    """

    df = stats[stats_name]
    
    min_y = df[stats_name].min()
    max_y = df[stats_name].max()

    fig, ax = plt.subplots(1, 1, figsize=(10, 4), dpi=300)

    if min_y == max_y:
        df2 = pd.DataFrame({"Value": [min_y - 1, min_y, min_y + 1], "Count": [0, len(df), 0]})
        sns.histplot(
            data=df, ax=ax, color='lightblue', shrink=0.1, discrete=True, kde=False, stat='count'
        )
        sns.lineplot(
            data=df2, x="Value", y="Count", ax=ax
        )
        ax.set_xticks([min_y - 1, min_y, min_y + 1])
    else:
        sns.histplot(
            data=df, x=stats_name, ax=ax, discrete=True, kde=True, stat='count',
        )
        
    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.ticklabel_format(axis='y', style='plain')

    return fig

def plot_nucleotides(stats, nucleotides, plot_type):

    """
    Plot a violin or boxen plot for nucleotide content.
    """
    df = stats['Per sequence nucleotide content']
    df = df.fillna(0)
    
    fig, ax = plt.subplots(figsize=(10, 4), dpi=300)
    
    if plot_type == 'violin':
        sns.violinplot(data=df, ax=ax, inner='quartile', cut=0, order=nucleotides)
    elif plot_type == 'boxen':
        sns.boxenplot(data=df, ax=ax, showfliers=False)
    else:
        logging.error(f"Unsupported plot type: {plot_type}. Supported types are 'violin' and 'boxen'.")
    
    ax.set_title('Nucleotide Content Distribution', fontsize=16)
    ax.set_ylabel('Frequency', fontsize=14)
    ax.set_xlabel('Nucleotides', fontsize=14)

    return fig

def plot_dinucleotides(stats, nucleotides, plot_type):
    """
    Plot a violin or boxen plot for dinucleotide content.
    """
    df = stats['Per sequence dinucleotide content']
    df = df.fillna(0)
    
    fig, axs = plt.subplots(len(nucleotides), 1, figsize=(10, len(nucleotides) * 3 + 2), sharey=True, dpi=300)

    for index, nt in enumerate(nucleotides):
        dinucleotides = [nt + nt2 for nt2 in nucleotides]
        df_nt = df[dinucleotides]
        if plot_type == 'violin':
            sns.violinplot(data=df_nt, ax=axs[index], inner='quartile', cut=0, order=dinucleotides)
        elif plot_type == 'boxen':
            sns.boxenplot(data=df_nt, ax=axs[index], showfliers=False, order=dinucleotides)
        else:
            logging.error(f"Unsupported plot type: {plot_type}. Supported types are 'violin' and 'boxen'.")
        
        if index == 0:
            axs[index].set_title('Dinucleotide content', fontsize=16)

        axs[index].set_xlabel('')
        axs[index].set_ylabel('Frequency', fontsize=14)
        axs[index].tick_params(axis='x', labelsize=12)
        axs[index].tick_params(axis='y', labelsize=12)

    axs[index].set_xlabel('Dinucleotide', fontsize=14)

    return fig

def plot_per_position_nucleotide_content(stats, stat_name, nucleotides, end_position=None):
    """
    Plot per position nucleotide content.
    Plot is constructed from 2 subplots:
    1. Nucleotide content per position.
    2. Number of sequences with length at least that position.
    The second plot is a line plot (with fill under it) showing how many sequences have length at least that position.
    If end_position is None, it will be set to the maximum sequence length.
    If end_position is specified, it will be used to limit the x-axis.
    """

    if not end_position:
        # Find the maximum sequence length from the stats
        seq_lengths = stats['Sequence lengths'].values.flatten()
        end_position = max(seq_lengths)
    else:
        # Ensure end_position is not greater than the maximum sequence length
        seq_lengths = stats['Sequence lengths'].values.flatten()
        max_length = int(max(seq_lengths))
        if end_position > max_length:
            logging.warning(f"end_position {end_position} is greater than the maximum sequence length {max_length}. Setting end_position to {max_length}.")
            end_position = max(seq_lengths)

    nucleotides = sorted(nucleotides)  # Ensure nucleotides are sorted for consistent plotting

    # Prepare the data for plotting
    # Create a DataFrame from the stats for the specified stat_name
    df = stats[stat_name]
    df = df.fillna(0)  # Fill NaN values with 0 for plotting

    # Ensure the DataFrame has the correct number of rows
    if end_position > len(df):
        # If end_position is greater than the number of rows in df, extend df with zeros
        additional_rows = pd.DataFrame(0, index=range(len(df), end_position), columns=df.columns)
        df = pd.concat([df, additional_rows], ignore_index=True)

    # Limit the DataFrame to the specified end_position
    df = df.iloc[:end_position]

    # Create the figure and axis for the plot
    fig, axs = plt.subplots(2, 1, figsize=(10, 4), height_ratios=[3, 1], dpi=300)

    # Plot the per position nucleotide content
    for nucleotide in nucleotides:
        axs[0].plot(df.index, df[nucleotide], label=nucleotide, linewidth=2)

    # set legend to the right side of the plot, outside the plot area
    axs[0].legend(title='Nucleotides', loc='upper right', fontsize=12, title_fontsize=14, bbox_to_anchor=(1.2, 1))
    axs[0].set_ylabel('Nucleotide content', fontsize=14)
    axs[0].set_title(f"{stat_name}", fontsize=16)

    # Plot the number of sequences with length at least that position
    length_counts = [sum(1 for length in seq_lengths if length >= pos) for pos in range(end_position)]
    # normalize length_counts to [0, 1]
    if length_counts:
        length_counts = [count / max(length_counts) for count in length_counts]

    axs[1].fill_between(range(end_position), length_counts, color='lightblue', alpha=0.5)
    axs[1].plot(range(end_position), length_counts, color='lightblue', linewidth=2)
    axs[1].set_xlabel('Position in sequence', fontsize=14)
    axs[1].set_ylabel('Proportion of\nsequences', fontsize=14)
    axs[1].yaxis.set_label_position("right")
    axs[1].yaxis.tick_right()
    axs[1].set_xticks(range(0, end_position + 1, max(1, end_position // 10)))
    axs[1].set_xticklabels(range(0, end_position + 1, max(1, end_position // 10)), fontsize=12)
    axs[1].set_ylim(0, max(length_counts) * 1.1)  # Set y-axis limits to [0, max(length_counts) * 1.1]  
    # set y-tick to 0, 0.5 and 1
    axs[1].set_yticks([0, 0.5, 1])
    axs[1].set_yticklabels(['0', '0.5', '1'], fontsize=12)

    # Adjust layout to remove extra space between subplots
    plt.subplots_adjust(wspace=0, hspace=0)

    return fig
