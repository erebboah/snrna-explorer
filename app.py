# app.py

import streamlit as st
from utils.routing import load_available_datasets, load_selected_metadata
from modules.metadata_explorer import metadata_explorer
from modules.gene_expression_explorer import gene_expression_explorer

st.markdown(
    "<h2 style='text-align: center;'>🧬 Welcome to the snRNA-seq Data Explorer</h2>",
    unsafe_allow_html=True,
)

st.markdown(
    "<p style='text-align: center;'>Select a dataset to begin exploring proportions and gene expression.</p>",
    unsafe_allow_html=True,
)

# Center the dropdown
col1, col2, col3 = st.columns([1, 2, 1])
with col2:
    dataset_options = ["Select a Dataset"] + load_available_datasets()
    dataset = st.selectbox("Choose a Dataset", dataset_options, label_visibility="collapsed")

if dataset == "Select a Dataset":
    st.stop()


# Load metadata once dataset is selected
metadata = load_selected_metadata(dataset)

# Horizontal divider
st.markdown("---")

# Section 1: Metadata proportions
st.subheader("📊 Metadata Proportion Explorer")
metadata_explorer(metadata)

# Horizontal divider
st.markdown("---")

# Section 2: Gene expression
st.subheader("🔥 Gene Expression Explorer")
gene_expression_explorer(dataset)
