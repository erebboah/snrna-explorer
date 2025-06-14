import streamlit as st
import pandas as pd
import os

def gene_expression_explorer(dataset):
    st.header("Gene Expression Explorer")
    st.markdown(f"Dataset selected: `{dataset}`")

    # Load gene reference table
    gene_df = pd.read_csv("data/genes.csv")
    gene_df = gene_df.dropna(subset=["gene_name", "gene_id"]).drop_duplicates()

    # Build searchable list
    all_genes = sorted(set(gene_df["gene_name"]) | set(gene_df["gene_id"]))

    # Input method: free text + optional autocomplete dropdown
    st.markdown("#### Select Genes of Interest")
    st.markdown("Paste gene names or IDs (comma- or space-separated). Autocomplete supported.")
    
    gene_input = st.text_input("Enter gene names/IDs", placeholder="e.g., Acta1, Myog, Dmd")

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

    # TODO: Add grouping options and data loading based on genes selected
