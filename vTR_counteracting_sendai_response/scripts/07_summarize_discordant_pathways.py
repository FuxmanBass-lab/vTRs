import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, leaves_list
import numpy as np

# === Load NES summary table ===
df = pd.read_csv("../results/gsea/summary_reactome_nes.tsv", sep="\t", index_col=0)

# Identify relevant columns
ref_col = "DE_all_stim_vector_vs_unstim_vector"
vtr_cols = [col for col in df.columns if re.match(r"DE_all_stim_VTR\d+_vs_vector", col)]

# Categorization function
def categorize(stim, vtr):
    if pd.isna(stim) or pd.isna(vtr):
        return None
    if vtr > 0 and stim < 0:
        return "discordant_up"
    elif vtr < 0 and stim > 0:
        return "discordant_down"
    else:
        return None

# Build discordant summary
summary = []
pathway_labels = {}

for pathway in df.index:
    stim_val = df.loc[pathway, ref_col]
    up_list, down_list = [], []
    for vtr_col in vtr_cols:
        vtr_val = df.loc[pathway, vtr_col]
        result = categorize(stim_val, vtr_val)
        vtr_id = re.search(r"VTR\d+", vtr_col).group(0)
        if result == "discordant_up":
            up_list.append(vtr_id)
        elif result == "discordant_down":
            down_list.append(vtr_id)
    total = len(up_list) + len(down_list)
    if total > 0:
        label = "discordant_up" if len(up_list) >= len(down_list) else "discordant_down"
        pathway_labels[pathway] = label
        summary.append({
            "Pathway": pathway,
            "Discordant_Up": len(up_list),
            "Discordant_Down": len(down_list),
            "Total_Discordant": total,
            "Up_VTRs": ", ".join(up_list),
            "Down_VTRs": ", ".join(down_list),
        })

# Save summary
discordant_df = pd.DataFrame(summary).sort_values(by="Total_Discordant", ascending=False)
discordant_df.to_csv("../results/gsea/discordant_pathways_summary.tsv", sep="\t", index=False)
print("Saved summary to ../results/gsea/discordant_pathways_summary.tsv")

# Build NES matrix
clean_vtr_ids = sorted(set(re.search(r"VTR\d+", col).group(0) for col in vtr_cols))
nes_matrix = pd.DataFrame(0.0, index=discordant_df["Pathway"], columns=clean_vtr_ids)

for _, row in discordant_df.iterrows():
    pathway = row["Pathway"]
    up_vtrs = row["Up_VTRs"].split(", ") if row["Up_VTRs"] else []
    down_vtrs = row["Down_VTRs"].split(", ") if row["Down_VTRs"] else []
    for vtr_id in up_vtrs + down_vtrs:
        vtr_col = f"DE_all_stim_{vtr_id}_vs_vector"
        if vtr_col in df.columns:
            nes_matrix.loc[pathway, vtr_id] = df.loc[pathway, vtr_col]

# Drop rows and columns with all 0
nes_matrix = nes_matrix.loc[~(nes_matrix == 0).all(axis=1), ~(nes_matrix == 0).all(axis=0)]

# === Clustermap ===
sns.set(font_scale=0.9, style='white')
cmap = sns.diverging_palette(240, 10, as_cmap=True)
vmin, vmax = -3, 3

g = sns.clustermap(
    nes_matrix,
    metric="cosine",
    method="average",
    cmap=cmap,
    vmin=vmin,
    vmax=vmax,
    figsize=(14, max(6, 0.3 * nes_matrix.shape[0])),
    dendrogram_ratio=(0.15, 0.03),
    cbar_pos=(0.98, 0.2, 0.02, 0.6),  # full-height colorbar, far right
    xticklabels=True,
    yticklabels=True
)

g.ax_heatmap.set_xlabel("VTRs")
g.ax_heatmap.set_ylabel("Pathways")
g.ax_heatmap.tick_params(axis='x', labelrotation=90)
g.fig.suptitle("Discordant Pathways NES by VTR", fontsize=14, y=1.03)

# Save
output_dir = "../results/"
os.makedirs(output_dir, exist_ok=True)
g.savefig(os.path.join(output_dir, "discordant_pathway_heatmap_clustermap.pdf"))
g.savefig(os.path.join(output_dir, "discordant_pathway_heatmap_clustermap.png"), dpi=300)

print("Saved final clustermap with aligned dendrogram and full-height colorbar.")