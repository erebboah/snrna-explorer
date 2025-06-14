import os
import json
import seaborn as sns
import streamlit as st

def get_color_dict(var_key, metadata):
    palette_path = f"color_palettes/{var_key}_palette.json"
    if os.path.exists(palette_path):
        with open(palette_path) as f:
            return json.load(f)
    return None

def show_color_picker(var_key, metadata, color_dict=None):
    st.markdown("#### Customize Colors for Groups")
    os.makedirs("color_palettes", exist_ok=True)
    unique_vals = sorted(metadata[var_key].dropna().unique())
    color_dict = color_dict or {}

    for i, val in enumerate(unique_vals):
        default = color_dict.get(val, "#%02x%02x%02x" % tuple(int(c * 255) for c in sns.color_palette("husl", len(unique_vals))[i]))
        picked = st.color_picker(f"{val}", default, key=f"color_{val}")
        color_dict[val] = picked

    if st.button("Save Color Palette"):
        with open(f"color_palettes/{var_key}_palette.json", "w") as f:
            json.dump(color_dict, f, indent=2)
        st.success(f"Palette saved to color_palettes/{var_key}_palette.json!")

    return color_dict

