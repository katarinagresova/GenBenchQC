import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def plot_nucleotides(stats1, stats2, result, dist_thresh, nucleotides, plot_type='violin'):
    """
    Plot the nucleotide content of two sequences.
    """
    if plot_type == 'violin':
        return violin_plot_nucleotides(stats1, stats2, result, dist_thresh, nucleotides)
    else:
        raise ValueError(f"Unknown plot type: {plot_type}")
    
def plot_dinucleotides(stats1, stats2, result, dist_thresh, nucleotides, plot_type='violin'):
    """
    Plot the dinucleotide content of two sequences.
    """
    if plot_type == 'violin':
        return violin_plot_dinucleotides(stats1, stats2, result, dist_thresh, nucleotides)
    else:
        raise ValueError(f"Unknown plot type: {plot_type}")

def plot_lengths(stats1, stats2, result, dist_thresh, plot_type='violin'):
    """
    Plot the sequence lengths of two sequences.
    """
    if plot_type == 'violin':
        return violin_plot_one_stat(stats1, stats2, 'Sequence lengths', result, dist_thresh, x_label='Sequence length', title='Sequence Length Distribution')
    else:
        raise ValueError(f"Unknown plot type: {plot_type}")
    
def plot_gc_content(stats1, stats2, result, dist_thresh, plot_type='violin'):
    """
    Plot the GC content of two sequences.
    """
    if plot_type == 'violin':
        return violin_plot_one_stat(stats1, stats2, 'Per sequence GC content', result, dist_thresh, x_label='GC content', title='GC Content Distribution')
    else:
        raise ValueError(f"Unknown plot type: {plot_type}")

def violin_plot_nucleotides(stats1, stats2, result, dist_thresh, nucleotides):

    df = melt_stats(stats1, stats2, 'Per sequence nucleotide content', var_name='Nucleotide', value_name='Frequency')
    min_y = df['Frequency'].min()
    max_y = df['Frequency'].max()

    fig, ax = plt.subplots(1, 1, figsize=(10, 4), dpi=300)
    sns.violinplot(
        x='Nucleotide', 
        y='Frequency', 
        hue="label", 
        split=True, 
        data=df[df['Nucleotide'].isin(nucleotides)],
        gap=.1,
        order=nucleotides,
        hue_order=[stats1.label, stats2.label],
        density_norm='width',
        palette=HuePalette(),
    )

    red_flag = False
    for index, nt in enumerate(nucleotides):
        if result[0][nt] > dist_thresh:
            red_flag = True
            # draw a red rectangle around the violins
            ax.add_patch(make_red_flag_rectangle(index, min_y, max_y))

    ax.set_title('Nucleotide content', fontsize=16)
    ax.set_xlabel('Nucleotide', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax = prepare_legend(ax, red_flag, dist_thresh)

    return fig


def violin_plot_dinucleotides(stats1, stats2, result, dist_thresh, nucleotides):

    df = melt_stats(stats1, stats2, 'Per sequence dinucleotide content', var_name='Dinucleotide', value_name='Frequency')
    min_y = df['Frequency'].min()
    max_y = df['Frequency'].max()
    
    fig, axs = plt.subplots(len(nucleotides), 1, figsize=(10, len(nucleotides) * 3 + 2), sharey=True, dpi=300)
    red_flag = False
    for index, nt in enumerate(nucleotides):
        dinucleotides = [nt + nt2 for nt2 in nucleotides]
        row = df[df['Dinucleotide'].isin(dinucleotides)]

        sns.violinplot(
            x='Dinucleotide', 
            y='Frequency', 
            hue="label", 
            split=True, 
            data=row,
            gap=.1,
            order=dinucleotides,
            hue_order=[stats1.label, stats2.label],
            ax=axs[index],
            density_norm='width',
            palette=HuePalette()
        )
        
        if index == 0:
            axs[index].set_title('Dinucleotide content', fontsize=16)

        axs[index].set_xlabel('')
        axs[index].legend().set_visible(False)
        axs[index].set_ylabel('Frequency', fontsize=14)
        axs[index].tick_params(axis='x', labelsize=12)
        axs[index].tick_params(axis='y', labelsize=12)

        for di_index, dint in enumerate(dinucleotides):
            if result[0][dint] > dist_thresh:
                red_flag = True
                # draw a red rectangle around the violins, put it behind the violins
                axs[index].add_patch(make_red_flag_rectangle(di_index, min_y, max_y))

    axs[index] = prepare_legend(axs[index], red_flag, dist_thresh)
    axs[index].set_xlabel('Dinucleotide', fontsize=14)

    return fig

def violin_plot_one_stat(stats1, stats2, stats_name, result, dist_thresh, x_label='', title=''):
    """
    Plot a single statistic from two stats objects.
    """

    # make dataframe with two columns: label and values
    df1 = pd.DataFrame(stats1.stats[stats_name].values(), columns=[stats_name])
    df2 = pd.DataFrame(stats2.stats[stats_name].values(), columns=[stats_name])
    df1['label'] = str(stats1.label)
    df2['label'] = str(stats2.label)
    df = pd.concat([df1, df2], ignore_index=True)
    
    min_y = df[stats_name].min()
    max_y = df[stats_name].max()

    fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=300)
    sns.violinplot(
        y=stats_name, 
        hue="label", 
        split=True, 
        data=df,
        gap=.1,
        hue_order=[str(stats1.label), str(stats2.label)],
        ax=ax,
        density_norm='width',
        palette=HuePalette()
    )
    
    # result is a tuple of (distances, passed)
    red_flag = False
    if result and result[0] > dist_thresh:
        red_flag = True
        # draw a red rectangle around the violins
        ax.add_patch(make_red_flag_rectangle(0, min_y, max_y))

    ax.set_title(title, fontsize=16)
    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel(stats_name, fontsize=14)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    if min_y == max_y:
        # show only one value
        ax.set_yticks([min_y])
    ax.ticklabel_format(axis='y', style='plain')
    ax = prepare_legend(ax, red_flag, dist_thresh) 

    return fig

def plot_per_base_sequence_comparison(stats1, stats2, stats_name, result, p_value_thresh, nucleotides, end_position, x_label='', title=''):

    df1 = pd.DataFrame(stats1.stats[stats_name]).T
    df2 = pd.DataFrame(stats2.stats[stats_name]).T

    fig, axs = plt.subplots(len(nucleotides), 1, figsize=(10, len(nucleotides) * 2 + 2), sharey=True, dpi=300)
    red_flag = False
    for index, nt in enumerate(nucleotides):

        end_position = min(end_position, len(result[0][nt]))

        df1_base = df1[nt][:end_position]
        df2_base = df2[nt][:end_position]
        axs[index].plot(df1.index[:end_position], df1_base, label=f"{stats1.label}", color=HuePalette()[0])
        axs[index].plot(df2.index[:end_position], df2_base, label=f"{stats2.label}", color=HuePalette()[1])

        axs[index].set_ylim(-0.1, 1.1)
        axs[index].set_ylabel('Frequency', fontsize=14)
        axs[index].legend().set_visible(False)
        axs[index].tick_params(axis='x', labelsize=12)
        axs[index].tick_params(axis='y', labelsize=12)

        # add text to the plot with the nucleotide name
        axs[index].text(0.9, 0.85, f'Nucleotide: {nt}', ha='center', va='bottom', fontsize=14, transform=axs[index].transAxes)

        if title and index == 0:
            axs[index].set_title(f'{title}', fontsize=16)


        if result:
            for i in range(end_position):
                
                p_value = result[0][nt][i]
                if p_value < p_value_thresh:
                    axs[index].axvspan(i-0.45, i+0.45, facecolor='red', alpha=0.2)
                    red_flag = True

    axs[index].set_xlabel(x_label, fontsize=14)
    axs[index].ticklabel_format(axis='both', style='plain')
    axs[index] = prepare_legend(axs[index], red_flag, p_value_thresh, box_to_anchor=(0.5, -0.3))

    return fig

def plot_plot_basic_descriptive_stats(stats1, stats2):

    # initialize the figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 3), dpi=300)

    # hide the axes
    ax.axis('off')

    unique_bases1 = sorted(stats1.stats["Unique bases"])
    unique_bases2 = sorted(stats2.stats["Unique bases"])
    # create a table with the basic descriptive statistics
    table_data = [
        ["Filename", stats1.filename, stats2.filename],
        ["Label", stats1.label, stats2.label],
        ["Number of sequences", stats1.stats["Number of sequences"], stats2.stats["Number of sequences"]],
        ["Unique sequences", stats1.stats["Number of sequences left after deduplication"], stats2.stats["Number of sequences left after deduplication"]],
        ["Number of bases", stats1.stats["Number of bases"], stats2.stats["Number of bases"]],
        ["Unique bases", ', '.join(unique_bases1), ', '.join(unique_bases2)],
        ["%GC content", f"{stats1.stats['%GC content']:.2%}", f"{stats2.stats['%GC content']:.2%}"],
        
    ]
    table = ax.table(cellText=table_data, colLabels=None, cellLoc='center', loc='center')
    # set first column to left align
    for cell in table.get_celld().values():
        if cell.get_text() and cell.get_text().get_text() in [row[0] for row in table_data]:
            cell.set_text_props(ha='left', va='center')
        else:
            cell.set_text_props(ha='center', va='center')
    table.auto_set_font_size(False)
    table.set_fontsize(14)
    table.scale(0.7, 1.7)  # scale the table to make it more readable
    table.auto_set_column_width([0, 1, 2])

    # set the title of the plot
    ax.set_title('Basic Descriptive Statistics', fontsize=16)
    ax.set_frame_on(False)  # remove the frame around the table
    ax.set_xticks([])  # remove the x ticks
    ax.set_yticks([])  # remove the y ticks
    ax.set_xticklabels([])  # remove the x tick labels
    ax.set_yticklabels([])  # remove the y tick labels
    ax.spines['top'].set_visible(False)  # remove the top spine
    ax.spines['right'].set_visible(False)  # remove the right spine
    ax.spines['left'].set_visible(False)  # remove the left spine
    ax.spines['bottom'].set_visible(False)  # remove the bottom spine
    ax.set_xlim(0, 1)  # set the x limits to [0, 1]
    ax.set_ylim(0, 1)  # set the y limits to [0, 1]

    return fig

def plot_duplicates(result):

    if result[0] is None:
        return None
    else:
        # initialize the figure
        fig, ax = plt.subplots(1, 1, figsize=(10, 4), dpi=300)

        # hide the axes
        ax.axis('off')

        # create a table with the duplicates max top 10 duplicates
        table_data = []
        for seq in list(result[0])[:10]:
            if len(seq) > 50:
                seq = seq[:50] + '...'
            table_data.append([seq])
        if len(result[0]) > 10:
            table_data.append(['... and more (saving to file)'])

        table = ax.table(cellText=table_data, colLabels=None, cellLoc='center', loc='center')

        table.auto_set_font_size(False)
        table.set_fontsize(14)
        table.scale(0.7, 1.7)  # scale the table to make it more readable
        table.auto_set_column_width([0, 1])

    # set the title of the plot
    ax.set_title('Duplicate Sequences', fontsize=16)
    ax.set_frame_on(False)  # remove the frame around the table
    ax.set_xticks([])  # remove the x ticks
    ax.set_yticks([])  # remove the y ticks
    ax.set_xticklabels([])  # remove the x tick labels
    ax.set_yticklabels([])  # remove the y tick labels
    ax.spines['top'].set_visible(False)  # remove the top spine
    ax.spines['right'].set_visible(False)  # remove the right spine
    ax.spines['left'].set_visible(False)  # remove the left spine
    ax.spines['bottom'].set_visible(False)  # remove the bottom spine
    ax.set_xlim(0, 1)  # set the x limits to [0, 1]
    ax.set_ylim(0, 1)  # set the y limits to [0, 1] 

    return fig


def melt_stats(stats1, stats2, stats_name, var_name='Metric', value_name='Value', keep_positions=False):
    """
    Melt the stats DataFrame to long format and add a label column.
    """
    df1 = pd.DataFrame(stats1.stats[stats_name]).T
    df1 = df1.melt(value_vars=df1.columns, var_name=var_name, value_name=value_name)
    df1['label'] = stats1.label

    df2 = pd.DataFrame(stats2.stats[stats_name]).T
    df2 = df2.melt(value_vars=df2.columns, var_name=var_name, value_name=value_name)
    df2['label'] = stats2.label

    df = pd.concat([df1, df2], ignore_index=True)

    return df

def prepare_legend(ax, red_flag, dist_thresh, box_to_anchor=(0.5, -0.2)):
    """
    Prepare the legend for the plot.
    """
    legend_handles = ax.get_legend_handles_labels()[0]
    legend_labels = ax.get_legend_handles_labels()[1]
    if red_flag:
        legend_handles += [plt.Rectangle((0, 0), 1, 1, color='red', alpha=0.2)]
        legend_labels += [f'Distance > {dist_thresh}']
    ax.legend(
        handles = legend_handles,
        labels = legend_labels,
        title='Label',
        title_fontsize='14',
        fontsize='12',
        loc='upper center', 
        bbox_to_anchor=box_to_anchor, 
        ncol=3,
    )
    return ax

def make_red_flag_rectangle(index, min_y, max_y, margin=0.02):
    """
    Create a red rectangle to highlight the violins that are above the threshold.
    """
    flag_box_width = 1 - 2 * margin
    flag_box_hight = max_y - min_y + 2 * margin
    return plt.Rectangle((index - 0.5 + margin, min_y - margin), flag_box_width, flag_box_hight, color='red', alpha=0.2, zorder=-1)

class HuePalette:

    _palette = None

    def __new__(palette):
        if palette._palette is None:
            palette._palette = sns.color_palette()[:2]

        return palette._palette
                
