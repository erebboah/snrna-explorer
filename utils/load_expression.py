import os
import scanpy as sc
import pandas as pd
import numpy as np


def load_and_extract_genes(dataset, metadata, selected_genes, groupby_col):
    data_dir = f"data/pseudobulk_by_genotype/{dataset}"
    needed_genotypes = sorted(metadata["Genotype"].dropna().unique())

    adata_list = []

    for genotype in needed_genotypes:
        adata_path = f"{data_dir}/{dataset}_{genotype}_pseudobulk.h5ad"
        if not os.path.exists(adata_path):
            print(f"Warning: File not found for {genotype}: {adata_path}")
            continue

        adata = sc.read_h5ad(adata_path)

        # Filter genes by gene_name or gene_id
        gene_mask = (
            adata.var.get("gene_name", pd.Series(index=adata.var.index)).isin(selected_genes) |
            adata.var.get("gene_id", pd.Series(index=adata.var.index)).isin(selected_genes)
        )

        if gene_mask.sum() == 0:
            print(f"No selected genes found in {genotype}")
            continue

        adata_sub = adata[:, gene_mask].copy()

        # Determine correct cluster key
        cluster_key = "leiden_R" if "leiden_R" in adata_sub.obs.columns else "leiden"

        # Filter adata_sub to matching lab_sample_id and cluster_key from metadata
        valid_sample_ids = set(metadata["lab_sample_id"].dropna().unique())
        valid_clusters = set(metadata[cluster_key].dropna().unique())

        adata_sub = adata_sub[adata_sub.obs["lab_sample_id"].isin(valid_sample_ids) & adata_sub.obs[cluster_key].isin(valid_clusters)].copy()

        if adata_sub.n_obs == 0:
            print(f"No matching observations found for {genotype} after filtering.")
            continue

        # Group and sum expression by groupby_col
        adata_sub.obs[groupby_col] = adata_sub.obs[groupby_col].astype(str)
        group_ids = adata_sub.obs[groupby_col].values
        unique_groups = np.unique(group_ids)

        X = []
        obs_names = []
        for group in unique_groups:
            mask = group_ids == group
            summed = adata_sub.X[mask].sum(axis=0)

            if hasattr(summed, "A1"):  # sparse matrix handling
                summed = summed.A1
            else:
                summed = np.asarray(summed).flatten()

            X.append(summed)
            obs_names.append(group)

        X = np.array(X)

        # Build new adata
        new_adata = sc.AnnData(X)
        new_adata.obs[groupby_col] = obs_names
        new_adata.obs['Genotype'] = genotype
        new_adata.var = adata_sub.var.loc[adata_sub.var_names]

        adata_list.append(new_adata)

    # Merge all
    if adata_list:
        merged = adata_list[0].concatenate(*adata_list[1:])
        #merged.write_h5ad('test.h5ad')
        return merged
    else:
        return None
