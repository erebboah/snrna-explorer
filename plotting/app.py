# app.py

import streamlit as st
import sqlite3
import pandas as pd
from plotting.stacked_barplot import stacked_barplot_proportions

st.set_page_config(layout="wide")
st.title("snRNA-seq Metadata Proportion Explorer")

@st.cache_data
def load_metadata():
    conn = sqlite3.connect("metadata/metadata.sqlite")
    df = pd.read_sql_query("SELECT * FROM metadata", conn)
    conn.close()
    return df

metadata = load_metadata()

# Show basic info
st.markdown("### Preview of Metadata")
st.dataframe(metadata.head())

# Select plotting keys
categorical_cols = metadata.select_dtypes(include="object").columns.tolist()
cluster_key = st.selectbox("Select Cluster Key (e.g., leiden)", categorical_cols)
var_key = st.selectbox("Select Grouping Variable (e.g., sex, genotype)", categorical_cols)

# Optional color dict
if st.checkbox("Use custom colors?"):
    unique_vals = sorted(metadata[var_key].dropna().unique())
    default_colors = sns.color_palette("husl", n_colors=len(unique_vals))
    color_dict = {val: default_colors[i] for i, val in enumerate(unique_vals)}
else:
    color_dict = None

# Plot button
if st.button("Generate Stacked Bar Plot"):
    fig = stacked_barplot_proportions(
        metadata, cluster_key=cluster_key, var_key=var_key,
        color_dict=color_dict
    )
    st.pyplot(fig)

