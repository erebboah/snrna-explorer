import streamlit as st
import pandas as pd
import scanpy as sc
import anndata as ad
import os
from utils.load_expression import load_and_extract_genes
from plotting.gene_heatmap import plot_gene_heatmap


def gene_expression_explorer(dataset, metadata):
    st.header("Gene Expression Heatmap")
    st.markdown(f"Dataset selected: `{dataset}`")

    # Load gene reference table
    gene_df = pd.read_csv("data/genes.csv")
    gene_df = gene_df.dropna(subset=["gene_name", "gene_id"]).drop_duplicates()

    # Build searchable list
    all_genes = sorted(set(gene_df["gene_name"]) | set(gene_df["gene_id"]))

    # Input method: free text + optional autocomplete dropdown
    st.markdown("### Select Genes of Interest")
    st.markdown("Paste gene names or IDs (comma- or space-separated).")
    
    gene_input = st.text_input("Enter gene names/IDs", placeholder="e.g., Slc17a7, Myog, Flt1")

    # Optional: multiselect dropdown for visual browsing
    selected_from_dropdown = st.multiselect(
        "Or select from list (optional)",
        options=all_genes,
        key="gene_dropdown"
    )

    # Process gene_input text field
    input_genes = [
        g.strip() for g in gene_input.replace(",", " ").split()
        if g.strip() in all_genes
    ]

    # Merge with dropdown-selected genes
    selected_genes = sorted(set(input_genes + selected_from_dropdown))

    if selected_genes:
        st.success(f"Selected {len(selected_genes)} gene(s): {', '.join(selected_genes)}")
    else:
        st.info("No valid genes selected yet.")


    # Determine cluster key
    cluster_key = "leiden_R" if "leiden_R" in metadata.columns else "leiden"

    # Identify columns that are constant within each sample
    sample_meta = (
        metadata.groupby("lab_sample_id", observed=True)
        .agg(lambda x: x.unique()[0] if len(set(x)) == 1 else pd.NA)
        .dropna(axis=1)
    )
    sample_cols = sample_meta.columns.tolist()

    # Identify columns that are constant within each cluster
    cluster_meta = (
        metadata.groupby(cluster_key, observed=True)
        .agg(lambda x: x.unique()[0] if len(set(x)) == 1 else pd.NA)
        .dropna(axis=1)
    )
    cluster_cols = cluster_meta.columns.tolist()

    # Union of valid grouping options, include keys used for aggregation
    # Desired keys if present
    must_include = ["lab_sample_id", "leiden", "leiden_R"]

    # Get valid columns from metadata
    existing_keys = [col for col in must_include if col in metadata.columns]

    # Combine with sample- and cluster-level grouping vars
    all_grouping_vars = existing_keys + [col for col in metadata.columns if col in set(sample_cols + cluster_cols) and col not in existing_keys]


    st.markdown("### Filter Data (Optional)")

    # Use this to track filtered metadata
    filtered_metadata = metadata.copy()

    # Allow user to add filters one by one from valid groupable columns
    with st.expander("Add filters"):
        filter_vars = st.multiselect("Select metadata fields to filter by:", all_grouping_vars)

        for filter_col in filter_vars:
            vals = sorted(filtered_metadata[filter_col].dropna().unique())
            selected = st.multiselect(f"Filter values for `{filter_col}`", vals, default=vals, key=f"filter_{filter_col}")
            filtered_metadata = filtered_metadata[filtered_metadata[filter_col].isin(selected)]

    st.markdown("### Choose groups")

    metadata = filtered_metadata

    # Step 1: choose grouping variable
    groupby_col = st.selectbox("Group samples by:", all_grouping_vars)

    # Step 2: infer unique group values from the filtered metadata
    unique_vals = sorted(metadata[groupby_col].dropna().unique())

    # Summary display only (no re-selection needed)
    st.markdown(f"**You selected {len(unique_vals)} group(s) under `{groupby_col}`**")


    # Confirm selected_genes is not empty before proceeding
    if selected_genes:
        merged_adata = load_and_extract_genes(dataset, metadata, selected_genes, groupby_col)

        if not merged_adata:
            st.warning("No matching gene(s) found in any available pseudobulk files.")
        else:
            loaded_genotypes = merged_adata.obs["Genotype"].unique().tolist()
            st.success(f"Loaded expression data for {len(loaded_genotypes)} genotype(s): {', '.join(loaded_genotypes)}")

            # Optional: preview underlying matrix
            df_preview = merged_adata.to_df().T
            df_preview.columns = merged_adata.obs[groupby_col]
            st.markdown("#### Expression Matrix Preview")
            st.dataframe(df_preview)

            # Plot heatmap if data was loaded
            fig = plot_gene_heatmap(merged_adata, groupby_col)
            st.pyplot(fig)
    else:
        st.info("Please enter one or more genes to proceed.")

