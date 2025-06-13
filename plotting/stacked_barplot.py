# plotting/stacked_barplot.py

import seaborn as sns
import matplotlib.pyplot as plt

def stacked_barplot_proportions(
    metadata, cluster_key, var_key, fsize=(5, 12), annotations=True,
    reverse_order=False, color_dict=None, cluster_order=None
):
    """
    Creates a horizontal stacked bar plot of proportions for a categorical variable across clusters.
    """
    grouped_data = metadata.groupby([cluster_key, var_key]).size().unstack(fill_value=0)
    proportions = grouped_data.div(grouped_data.sum(axis=1), axis=0)

    if cluster_order:
        proportions = proportions.loc[cluster_order]
    elif reverse_order:
        proportions = proportions.iloc[::-1]

    cluster_sizes = grouped_data.sum(axis=1)
    cluster_sizes = cluster_sizes.loc[proportions.index]

    unique_categories = grouped_data.columns.tolist()
    if color_dict:
        colors = [color_dict.get(cat, "gray") for cat in unique_categories]
    else:
        colors = sns.color_palette("husl", n_colors=len(unique_categories))

    fig, ax = plt.subplots(figsize=fsize)
    proportions.plot(kind="barh", stacked=True, color=colors, ax=ax, width=0.8, edgecolor=None)

    if annotations:
        for i, txt in enumerate(cluster_sizes):
            ax.text(1.02, i, f"{txt:,}", fontsize=12, va="center", transform=ax.get_yaxis_transform())

    ax.set_xlim(0, 1.1)
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    ax.set_xlabel("Proportion", fontsize=14)
    ax.set_ylabel(cluster_key, fontsize=14)
    ax.set_title(f'{var_key} by {cluster_key}', fontsize=16)

    if annotations:
        ax.legend(title=var_key, bbox_to_anchor=(1.2, 1), loc="upper left")
    else:
        ax.get_legend().remove()

    plt.grid(False)
    return fig

