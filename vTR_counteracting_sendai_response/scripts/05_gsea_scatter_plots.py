#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import re

# === Load NES summary ===
df = pd.read_csv("../results/gsea/summary_reactome_nes.tsv", sep="\t")

# Reference (stim vs unstim)
ref_col = "DE_all_stim_vector_vs_unstim_vector"
if ref_col not in df.columns:
    raise ValueError(f"Missing column '{ref_col}'")

# All VTRs
vtr_cols = [col for col in df.columns if re.match(r"DE_all_stim_VTR\d+_vs_vector", col)]
print(f"Found {len(vtr_cols)} VTRs")

# Output file
output_path = "../results/gsea/vtr_scatter_nes.pdf"
os.makedirs(os.path.dirname(output_path), exist_ok=True)

# === Category logic and color mapping ===
category_map = {
    0: ("#C0C0C0", "Stim-only regulated"),
    1: ("#8B1C3B", "VTR↑ & Stim↓"),
    2: ("#2166AC", "VTR↓ & Stim↑"),
    3: ("#8D8DFF", "Down in both"),
    4: ("#FF8080", "Up in both"),
    5: ("#FDE0DD", "VTR↑ only"),
    6: ("#D0E1F2", "VTR↓ only")
}

def categorize(stim, vtr):
    if pd.isna(stim) or pd.isna(vtr):
        return None
    if vtr == 0 and stim != 0:
        return 0  # Stim-only
    elif vtr > 0 and stim < 0:
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
        return None
    
    
# === Plot ===
with PdfPages(output_path) as pdf:
    saved = 0
    for vtr in vtr_cols:
        subset = df[[ref_col, vtr]].copy()
        subset.columns = ["Stim", "VTR"]
        subset = subset.dropna()


        subset["Stim_plot"] = subset["Stim"]
        subset["VTR_plot"] = subset["VTR"]

        subset["Category"] = [categorize(stim, vtr) for stim, vtr in zip(subset["Stim"], subset["VTR"])]
        subset = subset.dropna(subset=["Category"])

        if subset.empty:
            continue

        plt.figure(figsize=(5, 5))
        for cat_id, (color, label) in category_map.items():
            points = subset[subset["Category"] == cat_id]
            if not points.empty:
                plt.scatter(points["Stim_plot"], points["VTR_plot"],
                            c=color, label=f"{label} (n={len(points)})",
                            s=20, edgecolors='black', linewidths=0.3)

        plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
        plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)
        # plt.title(f"{vtr}", fontsize=9)
        plt.xlabel("Stimulated vs Unstimulated Vector (NES)", fontsize=8)
        vtr_id = re.search(r"VTR\d+", vtr).group()  # Extract "VTR46"
        plt.ylabel(f"Stimulated {vtr_id} vs Stimulated Vector (NES)", fontsize=8)

        # Place legend dynamically
        plt.legend(loc="best", fontsize=6, frameon=True)
        plt.tight_layout()
        pdf.savefig()
        plt.close()
        saved += 1

print(f"Done. Saved {saved} scatter plots to:\n{output_path}")