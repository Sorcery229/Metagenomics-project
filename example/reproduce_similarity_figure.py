"""Reproduce the similarity-distribution figure from committed data.

Mirrors the binning in notebooks/visualization.ipynb (cell reading
results/output_file.csv): per-pair similarity is n_of_1 / subject_length,
binned into categories, and the category counts are plotted as a bar chart.

Run from the repository root:
    python example/reproduce_similarity_figure.py

Outputs (written next to this script):
    example/similarity_distribution.png        - the bar chart
    example/similarity_counts.csv              - category -> count

Verify reproducibility by diffing similarity_counts.csv against the committed
example/expected_similarity_counts.csv (the deterministic part; the PNG may
differ at the pixel level across matplotlib versions).
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

df = pd.read_csv(os.path.join(ROOT, "results/output_file.csv"), dtype="str")
df["n_of_1"] = pd.to_numeric(df["n_of_1"], errors="coerce")
df["subject_length"] = pd.to_numeric(df["subject_length"], errors="coerce")
df["ratio"] = df["n_of_1"] / df["subject_length"]

def categorize_ratio(r):
    if r < 0.1:
        return "Below 10%"
    elif r < 0.7:
        return "Below 70%"
    elif r == 1.0:
        return "100%"
    else:
        return f"{int(r * 100)}%"

df["ratio_category"] = df["ratio"].apply(categorize_ratio)
counts = df["ratio_category"].value_counts().sort_index()
counts.rename_axis("category").rename("count").to_csv(
    os.path.join(HERE, "similarity_counts.csv"))

fig, ax = plt.subplots(figsize=(9, 4))
palette = sns.color_palette("RdYlGn_r", n_colors=len(counts))
ax.bar(counts.index.astype(str), counts.values, color=palette)
ax.set_xlabel("per-base similarity category")
ax.set_ylabel("number of genome pairs")
ax.set_title("Cross-database per-base similarity distribution")
for lbl in ax.get_xticklabels():
    lbl.set_rotation(60)
    lbl.set_ha("right")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "similarity_distribution.png"), dpi=150)
print("categories:", len(counts), "| total pairs:", int(counts.sum()))
print(counts.to_string())
