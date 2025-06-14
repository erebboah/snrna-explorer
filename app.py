# app.py

import streamlit as st
import os
from utils.metadata_loader import load_metadata
from utils.color_palette import get_color_dict, show_color_picker
from plotting.stacked_barplot import stacked_barplot_proportions

st.set_page_config(layout="wide")
st.title("snRNA-seq Data Explorer")

# Set metadata file path (can update later to be interactive)
data_path = "data/metadata/IGVFB01_LeftCortex_metadata.sqlite"

# Check if metadata file exists
if not os.path.exists(data_path):
    st.error(f"Metadata file not found at: {data_path}")
    st.stop()

# Load metadata
metadata = load_metadata(data_path)

# Show basic info
st.markdown("### Preview of Metadata")
st.dataframe(metadata.head())

# Get valid grouping columns (exclude individual cell IDs)
categorical_cols = metadata.select_dtypes(include="object").columns.tolist()
excluded_cols = ["cellID", "cell_id", "CellID"]
valid_group_cols = [
    col for col in categorical_cols
    if col not in excluded_cols and metadata[col].nunique() < len(metadata)
]

# User selects keys for plotting
cluster_key = st.selectbox("Select Cluster Key (e.g., leiden)", valid_group_cols)
var_key = st.selectbox("Select Grouping Variable (e.g., sex, genotype)", valid_group_cols)

# Load or define color palette
color_dict = get_color_dict(var_key, metadata)

# Optional color customization
if st.checkbox("Customize colors?"):
    color_dict = show_color_picker(var_key, metadata, color_dict)

# Plot button
if st.button("Generate Stacked Bar Plot"):
    fig = stacked_barplot_proportions(
        metadata,
        cluster_key=cluster_key,
        var_key=var_key,
        color_dict=color_dict
    )
    st.pyplot(fig)