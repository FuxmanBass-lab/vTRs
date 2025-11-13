#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import re
import os
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

# === Load NES summary table ===
df = pd.read_csv("../results/gsea/summary_reactome_nes.tsv", sep="\t", index_col=0)

# === Identify columns ===
ref_col = "DE_all_stim_vector_vs_unstim_vector"
vtr_cols = [col for col in df.columns if re.match(r"DE_all_stim_VTR\d+_vs_vector", col)]
print(f"Found {len(vtr_cols)} VTRs")

# === Categorization function ===
def categorize(stim, vtr):
    if pd.isna(stim) or pd.isna(vtr):
        return np.nan
    if vtr == 0:
        return 0
    if vtr > 0 and stim < 0:
        return 1
    elif vtr < 0 and stim > 0:
        return 2
    elif vtr < 0 and stim < 0:
        return 3
    elif vtr > 0 and stim > 0:
        return 4
    elif vtr > 0 and stim == 0:
        return 5
    elif vtr < 0 and stim == 0:
        return 6
    else:
        return np.nan

# === Build categorization matrix ===
cat_matrix = pd.DataFrame(index=df.index, columns=vtr_cols)
for col in vtr_cols:
    cat_matrix[col] = [
        categorize(df.loc[pathway, ref_col], df.loc[pathway, col])
        for pathway in df.index
    ]
cat_matrix = cat_matrix.astype(float).fillna(0)

# === Remove VTRs with no regulated pathways (all zeros) ===
nonzero_mask = (cat_matrix != 0).any(axis=0)
vtr_nonzero_cols = cat_matrix.columns[nonzero_mask]
print(f"{len(vtr_nonzero_cols)} VTRs kept after removing unregulated ones.")

# Clean up column names before selecting non-zero ones
cleaned_cols = {col: re.search(r"VTR\d+", col).group(0) for col in vtr_nonzero_cols}
cat_matrix = cat_matrix[vtr_nonzero_cols]
cat_matrix.rename(columns=cleaned_cols, inplace=True)

# === Compute CSI-based clustering ===
print("Clustering VTRs using Connection Specificity Index (CSI)...")
corr_matrix = cat_matrix.corr(method="pearson")
ranks = corr_matrix.rank(axis=1, ascending=False)
n = len(vtr_nonzero_cols)
csi_matrix = 1 - ranks / n
csi_distance = 1 - csi_matrix
distance_matrix = squareform(csi_distance.values, checks=False)
col_linkage = linkage(distance_matrix, method="average")

# === Load desired pathway order ===
with open("../data/pathways.txt") as f:
    ordered_pathways = [line.strip() for line in f if line.strip() in cat_matrix.index]

# Reorder rows based on the file
cat_matrix = cat_matrix.loc[ordered_pathways]

# === Define color palette and legend ===
category_colors = {
    0: "#FFFFFF",  # Not regulated
    1: "#8B1C3B",  # VTR↑ & Stim↓ (deep red)
    2: "#2166AC",  # VTR↓ & Stim↑ (deep blue)
    3: "#8D8DFF",  # Down in both
    4: "#FF8080",  # Up in both
    5: "#FDE0DD",  # VTR↑ only
    6: "#D0E1F2",  # VTR↓ only
}

cmap = ListedColormap([category_colors[i] for i in sorted(category_colors.keys())])
legend_labels = [
    Patch(facecolor=category_colors[0], edgecolor='black', label="0: Not Regulated in VTR"),
    Patch(facecolor=category_colors[1], edgecolor='black', label="1: VTR↑ & Stim↓"),
    Patch(facecolor=category_colors[2], edgecolor='black', label="2: VTR↓ & Stim↑"),
    Patch(facecolor=category_colors[3], edgecolor='black', label="3: Down in both"),
    Patch(facecolor=category_colors[4], edgecolor='black', label="4: Up in both"),
    Patch(facecolor=category_colors[5], edgecolor='black', label="5: VTR↑ only"),
    Patch(facecolor=category_colors[6], edgecolor='black', label="6: VTR↓ only"),
]

# === Plot heatmap ===
sns.set(font_scale=0.9, style='white')
g = sns.clustermap(
    cat_matrix,
    col_linkage=col_linkage,
    col_cluster=True,
    row_cluster=False,
    cmap=cmap,
    xticklabels=True,
    yticklabels=False,
    linewidths=0.1,
    linecolor='white',
    figsize=(14, 14),
    dendrogram_ratio=(.1, .1),
    cbar_pos=None,
)

g.ax_heatmap.xaxis.set_ticks_position('top')
g.ax_heatmap.xaxis.set_label_position('top')
g.ax_heatmap.tick_params(axis='x', labelrotation=90, labelsize=8)
g.fig.suptitle("VTR Pathway Regulation Patterns (CSI Clustered)", fontsize=14, y=1.04)

# === Add legend ===
g.ax_heatmap.legend(
    handles=legend_labels,
    title='Regulation Type',
    loc='lower center',
    bbox_to_anchor=(0.5, -0.3),
    ncol=3,
    frameon=False,
    fontsize=9,
    title_fontsize=10,
)

g.ax_heatmap.collections[0].set_rasterized(True)

# === Save ===
os.makedirs("../results", exist_ok=True)
g.fig.tight_layout()
g.fig.savefig("../results/vtr_pathway_heatmap_clustered_CSI.png", dpi=800)
g.fig.savefig("../results/vtr_pathway_heatmap_clustered_CSI.pdf", format="pdf", dpi=800)

# === Export matrix ===
col_order = g.dendrogram_col.reordered_ind
cat_matrix.iloc[:, col_order].to_csv("../results/vtr_pathway_heatmap_clustered_CSI.tsv", sep="\t")

print("CSI-clustered heatmap and matrix saved.")
