# utils/routing.py

import os
import sqlite3
import pandas as pd
import streamlit as st

def load_available_datasets():
    files = os.listdir("data/metadata/")
    
    # Extract base names (before "_metadata") for parquet and sqlite
    parquet_datasets = {
        f.replace("_metadata.parquet", "") for f in files if f.endswith(".parquet")
    }
    sqlite_datasets = {
        f.replace("_metadata.sqlite", "") for f in files if f.endswith(".sqlite")
    }

    # Use parquet if available, otherwise fallback to sqlite
    all_datasets = sorted(parquet_datasets.union(sqlite_datasets))
    return all_datasets

@st.cache_data
def load_selected_metadata(dataset):
    parquet_path = f"data/metadata/{dataset}_metadata.parquet"
    if os.path.exists(parquet_path):
        return pd.read_parquet(parquet_path)
    
    # fallback to sqlite
    db_path = f"data/metadata/{dataset}_metadata.sqlite"
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT * FROM metadata", conn)
    conn.close()
    return df

