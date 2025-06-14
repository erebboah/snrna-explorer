import streamlit as st
import pandas as pd
from utils.color_palette import get_color_dict, show_color_picker
from plotting.stacked_barplot import stacked_barplot_proportions

def metadata_explorer(metadata):
    st.markdown("### Preview of Metadata")
    st.dataframe(metadata.head())

    # Get valid grouping columns
    categorical_cols = metadata.select_dtypes(include="object").columns.tolist()
    excluded_cols = ["cellID", "cell_id", "CellID"]
    valid_group_cols = [
        col for col in categorical_cols
        if col not in excluded_cols and metadata[col].nunique() < len(metadata)
    ]

    # User selects cluster key and variable to group by
    cluster_key = st.selectbox("Select Cluster Key (e.g., leiden)", valid_group_cols)
    var_key = st.selectbox("Select Grouping Variable (e.g., sex, genotype)", valid_group_cols)

    # Load or generate color palette
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
