import sqlite3
import pandas as pd

def load_metadata(path):
    conn = sqlite3.connect(path)
    df = pd.read_sql_query("SELECT * FROM metadata", conn)
    conn.close()
    return df