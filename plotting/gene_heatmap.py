# plotting/gene_heatmap.py

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_gene_heatmap(adata, groupby_col):
    # Convert AnnData to DataFrame (genes x groups)
    df = adata.to_df().T
    df.columns = adata.obs[groupby_col]
    
    # Row-scale the data (genes)
    df = df.div(df.max(axis=1), axis=0)

    # Plot
    fig, ax = plt.subplots(figsize=(max(6, len(df.columns) * 0.6), max(4, len(df) * 0.4)))
    sns.heatmap(df, cmap="viridis", annot=False, linewidths=0.2, linecolor='gray', ax=ax)

    ax.set_xlabel("Sample Group")
    ax.set_ylabel("Gene")
    ax.set_title("Normalized Gene Expression Heatmap")

    plt.tight_layout()
    return fig
