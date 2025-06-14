import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import anndata as ad
import os
import sys
import re
import time
import sqlite3

if len(sys.argv) != 2:
    print("Usage: python make_pseudobulk_by_genotype.py <path_to_h5ad>")
    sys.exit(1)

input_path = sys.argv[1]
if not os.path.exists(input_path):
    print(f"File not found: {input_path}")
    sys.exit(1)

# Extract identifier from input file (strip extension + common suffix)
basename = os.path.basename(input_path)
identifier = re.sub(r'_annotated(_adata)?$', '', os.path.basename(input_path).replace('.h5ad', ''))
print(f"Processing: {identifier}")

start_time = time.time()

adata = sc.read_h5ad(input_path, backed = 'r')

# Choose cluster key
cluster_key = "leiden_R" if "leiden_R" in adata.obs.columns else "leiden"
if "lab_sample_id" not in adata.obs.columns or cluster_key not in adata.obs.columns:
    raise ValueError("adata.obs must contain both 'lab_sample_id' and a cluster key ('leiden' or 'leiden_R')")

# Group cell indices by lab_sample_id + cluster
grouped = adata.obs.groupby(["lab_sample_id", cluster_key], observed=True).indices

# Prepare summed count matrix and new .obs
X_rows = []
obs_records = []

for (sample_id, cluster_id), idx in grouped.items():
    counts = adata.layers["cellbender_counts"][idx].sum(axis=0)
    counts = counts.A1 if sparse.issparse(counts) else np.ravel(counts)
    X_rows.append(counts)
    obs_records.append({"lab_sample_id": sample_id, cluster_key: cluster_id})

X = sparse.csr_matrix(np.vstack(X_rows))
obs = pd.DataFrame(obs_records)
obs.index = [f"{row['lab_sample_id']}_{row[cluster_key]}" for _, row in obs.iterrows()]
var = adata.var.copy()

# Build pseudobulk AnnData
adata_pb = ad.AnnData(X=X, obs=obs, var=var)

# Collect sample-level metadata
sample_meta = (
    adata.obs
    .groupby("lab_sample_id", observed=True)
    .agg(lambda x: x.unique()[0] if len(set(x)) == 1 else np.nan)
    .dropna(axis=1)
)

# Collect cluster-level metadata
cluster_meta = (
    adata.obs
    .groupby(cluster_key, observed=True)
    .agg(lambda x: x.unique()[0] if len(set(x)) == 1 else np.nan)
    .dropna(axis=1)
)

# Determine overlapping column names
overlap_cols = sample_meta.columns.intersection(cluster_meta.columns)
cluster_meta_clean = cluster_meta.drop(columns=overlap_cols)

# Merge metadata
adata_pb.obs = adata_pb.obs.join(sample_meta, on="lab_sample_id")
adata_pb.obs = adata_pb.obs.join(cluster_meta_clean, on=cluster_key)

# Output directory
output_dir = os.path.join("pseudobulk_by_genotype", identifier)
os.makedirs(output_dir, exist_ok=True)

# Save one file per genotype
for genotype in adata_pb.obs["Genotype"].unique():
    adata_g = adata_pb[adata_pb.obs["Genotype"] == genotype].copy()
    out_path = os.path.join(output_dir, f"{identifier}_{genotype}_pseudobulk.h5ad")
    adata_g.write_h5ad(out_path, compression="gzip")
    print(f"Saved: {out_path}")


# Save original metadata to SQLite
os.makedirs("metadata", exist_ok=True)
sqlite_path = os.path.join("metadata", f"{identifier}_metadata.sqlite")

conn = sqlite3.connect(sqlite_path)
adata.obs.reset_index().to_sql("metadata", conn, if_exists="replace", index=False)
conn.close()

print(f"Saved metadata to {sqlite_path}")

elapsed = time.time() - start_time
print(f"Elapsed time: {elapsed:.2f} seconds.")