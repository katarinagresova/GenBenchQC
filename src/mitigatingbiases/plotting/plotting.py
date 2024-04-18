import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd

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

def plot_per_base_sequence_comparison(pos_df, neg_df, end_position, title='', x_label='', ax=None, stats=None, p_value_thresh=0.01):

    # get columns names
    bases = pos_df.columns.values
    # remove sum from bases (it doesn't have to be the last value)
    bases = bases[bases != 'sum']

    if ax is None:
        fig, ax = plt.subplots(nrows=len(bases), ncols=1, figsize=(15, 3 * len(bases)), sharex=True, sharey=True)

    for index, base in enumerate(bases):

        pos_freq = pos_df[base][:end_position] / pos_df['sum'][:end_position]
        neg_freq = neg_df[base][:end_position] / neg_df['sum'][:end_position]

        ax[index].plot(pos_df.index[:end_position], pos_freq, label=f"positive {base}")
        ax[index].plot(neg_df.index[:end_position], neg_freq, label=f"negative {base}")

        ax[index].set_ylim(-0.1, 1.1)
        ax[index].set_ylabel('Frequency')
        ax[index].legend()
        ax[index].set_title(f'Base: {base}')

        if stats:
            for i in range(len(pos_freq)):
                
                p_value = stats[base][i]
                if p_value < p_value_thresh:
                    ax[index].plot(i, p_value, 'ro',)

                    # plot salmon background on position with p-value < 0.05
                    ax[index].axvspan(i-0.45, i+0.45, facecolor='salmon', alpha=0.5)

                # add new value to legend
                ax[index].legend([f"positive {base}", f"negative {base}", f"p-value < 0.05"])

    ax[index].set_xlabel(x_label)

    return ax

def plot_composition_comparison(pos_df, neg_df, title='', x_label='', ax=None, stats=None, p_value_thresh=0.01):

    # get bcolumns names
    bases = pos_df.columns.values
    # remove sum from bases (it doesn't have to be the last value)
    bases = bases[bases != 'sum']


    if ax is None:
        fig, ax = plt.subplots(nrows=len(bases), ncols=1, figsize=(10, 2*len(bases)))

    #stats_passed = True
    for index, base in enumerate(bases):

        pos_freq = pos_df[base] / pos_df['sum']
        neg_freq = neg_df[base] / neg_df['sum']

        sns.histplot(pos_freq, ax=ax[index], label='positive', kde=True, alpha=0.3)
        sns.histplot(neg_freq, ax=ax[index], label='negative', kde=True, alpha=0.3)

        if stats:
            p_value = stats.get(base, 0)
            if p_value < p_value_thresh:
                ax[index].text(0.9, max( pos_df['sum'])/10, f"p-value: {p_value:.2e}", ha='center')
                
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
        ax[index].set_ylabel('Count')
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

def plot_lenght_comparison(pos_df, neg_df, title='', x_label='', ax=None, stats=None, p_value_thresh=0.01):

    if ax is None:
        fig, ax = plt.subplots(1, ncols=1, figsize=(10, 2))

    base = 'sum'

    pos_freq = pos_df[base]
    neg_freq = neg_df[base]

    sns.histplot(pos_freq, ax=ax, label='positive', kde=True, alpha=0.2, bins=20)
    sns.histplot(neg_freq, ax=ax, label='negative', kde=True, alpha=0.2, bins=20)

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
    ax.set_ylabel('Count')
    ax.legend()

    ax.set_xlabel(x_label)

    return ax

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