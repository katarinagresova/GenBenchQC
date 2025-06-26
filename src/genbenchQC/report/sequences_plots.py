import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def hist_plot_one_stat(stats, stats_name, x_label='', title=''):
    """
    Plot a single statistic from two stats objects.
    """

    df = pd.DataFrame(stats[stats_name].values(), columns=[stats_name])
    
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
        

    ax.set_xlabel(x_label)
    ax.set_ylabel('Count')
    ax.set_title(title)

    return fig
        