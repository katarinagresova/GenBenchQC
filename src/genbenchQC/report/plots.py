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


def violin_plot_nucleotides(stats1, stats2, result, dist_thresh, nucleotides):

    df = melt_stats(stats1, stats2, 'Per sequence nucleotide content', var_name='Nucleotide', value_name='Percentage')
    min_y = df['Percentage'].min()
    max_y = df['Percentage'].max()

    fig, ax = plt.subplots(1, 1, figsize=(10, 4), dpi=300)
    sns.violinplot(
        x='Nucleotide', 
        y='Percentage', 
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
    ax.set_ylabel('Percentage', fontsize=14)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax = prepare_legend(ax, red_flag, dist_thresh)

    return fig


def violin_plot_dinucleotides(stats1, stats2, result, dist_thresh, nucleotides):

    df = melt_stats(stats1, stats2, 'Per sequence dinucleotide content', var_name='Dinucleotide', value_name='Percentage')
    min_y = df['Percentage'].min()
    max_y = df['Percentage'].max()
    
    fig, axs = plt.subplots(len(nucleotides), 1, figsize=(10, len(nucleotides) * 3 + 2), sharey=True, dpi=300)
    red_flag = False
    for index, nt in enumerate(nucleotides):
        dinucleotides = [nt + nt2 for nt2 in nucleotides]
        row = df[df['Dinucleotide'].isin(dinucleotides)]

        sns.violinplot(
            x='Dinucleotide', 
            y='Percentage', 
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
        axs[index].set_ylabel('Percentage', fontsize=14)
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

def plot_per_base_sequence_comparison(stats1, stats2, stats_name, result, p_value_thresh, nucleotides, end_position, x_label='', title=''):

    df1 = pd.DataFrame(stats1.stats[stats_name]).T
    df2 = pd.DataFrame(stats2.stats[stats_name]).T

    fig, axs = plt.subplots(len(nucleotides), 1, figsize=(10, len(nucleotides) * 3 + 2), sharey=True, dpi=300)
    red_flag = False
    for index, nt in enumerate(nucleotides):

        end_position = min(end_position, len(df1[nt]), len(df2[nt]))

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
    axs[index] = prepare_legend(axs[index], red_flag, p_value_thresh)

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

def prepare_legend(ax, red_flag, dist_thresh):
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
        bbox_to_anchor=(0.5, -0.2), 
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
                
