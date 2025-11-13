#!/usr/bin/env python3
"""
generate_heatmap.py

Reads a tab-separated file of protein proportions and an annotation file,
aggregates by annotation group, orders the groups in a custom sequence,
and draws a heatmap.

Usage:
    python generate_heatmap.py data.tsv annotation.tsv output_heatmap.png
    python generate_heatmap_group.py ../data/ppi_proportions.tsv ../data/annots.tsv ../results/group_heatmap.png
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_style('white')
sns.set_context('notebook', font_scale=0.8)

def generate_heatmap(proportions_tsv, annotation_tsv, output_png):
    # 1. Load per‐protein proportions
    df = pd.read_csv(proportions_tsv, sep='\t', index_col='Protein')

    # 2. Load annotation (Protein → Group)
    annot = pd.read_csv(annotation_tsv, sep='\t', header=None,
                        names=['Protein','Group'], index_col='Protein')

    # 3. Merge and drop any proteins lacking an annotation
    df = df.join(annot, how='inner')

    # 4. Aggregate by Group (sum across proteins)
    group_df = df.groupby('Group').sum()
    # add data-driven pseudocount using half the smallest non-zero value
    values = group_df.values
    nonzeros = values[values > 0]
    min_nonzero = np.min(nonzeros) if nonzeros.size > 0 else 2e-1
    pseudocount = min_nonzero / 2.0
    group_df = group_df + pseudocount

    # 5a. Compute background proportion (mean across all conditions per group)
    background = group_df.mean(axis=1)
    # 5b. Compute fold-enrichment
    enrichment_df = group_df.div(background, axis=0)

    # 5. Reorder to your custom list
    desired_order = [
        'TF',
        'Coactivator',
        'Corepressor',
        'Other Cofactor',
        'RNA Processing',
        'Proteostasis',
        'Signaling',
        'Metabolism',
        'Immunity',
        'Other'
    ]
    # keep only those present, in that order
    present = [g for g in desired_order if g in enrichment_df.index]
    group_df = enrichment_df.loc[present]

    # convert to log2 fold-enrichment for diverging heatmap
    group_df = np.log2(group_df)
    # clip values to ±1.5 for visualization
    group_df = group_df.clip(lower=-1.5, upper=1.5)
    # fixed color range ±1.5
    vmin, vmax = -1.5, 1.5

    # 6. Plot heatmap
    plt.figure(figsize=(8, max(4, 0.4 * len(group_df))))
    ax = sns.heatmap(
        group_df,
        cmap='bwr',
        cbar_kws={'label': 'Log2 Fold Enrichment'},
        linewidths=0.2,
        linecolor='white',
        xticklabels=True,
        yticklabels=True,
        center=0,
        vmin=vmin,
        vmax=vmax
    )

    ax.set_xlabel('Condition')
    ax.set_ylabel('Annotation Group')
    ax.set_title('Protein Proportions by Annotation Group')
    plt.xticks(rotation=45, ha='right', fontsize='small')
    plt.tight_layout()

    # 7. Save
    plt.savefig(output_png, dpi=300)
    plt.close()

def main():
    if len(sys.argv) != 4:
        print("Usage: python generate_heatmap.py <data.tsv> <annotation.tsv> <output.png>")
        sys.exit(1)
    generate_heatmap(sys.argv[1], sys.argv[2], sys.argv[3])

if __name__ == '__main__':
    main()