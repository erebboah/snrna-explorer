# utils/routing.py

import os
import sqlite3
import pandas as pd
import streamlit as st

def load_available_datasets():
    files = os.listdir("data/metadata/")
    return sorted([f.replace("_metadata.sqlite", "") for f in files if f.endswith(".sqlite")])

@st.cache_data
def load_selected_metadata(dataset):
    db_path = f"data/metadata/{dataset}_metadata.sqlite"
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT * FROM metadata", conn)
    conn.close()
    return df
