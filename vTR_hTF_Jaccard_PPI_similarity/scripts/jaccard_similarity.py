#!/usr/bin/env python3
"""
Compute Jaccard similarity on a feature matrix.

- Input: a delimited table where rows are items (VTRs or HTFs) and columns are features (e.g., proteins). If you want column-wise similarities, pass --axis columns.
- Output: a square similarity matrix (CSV) with Jaccard similarities.

Examples
--------
# Auto-detect delimiter (comma/TSV), use first column as index, compute row-wise Jaccard
python jaccard_similarity.py matrix.csv --out jaccard_rows.csv

# Column-wise Jaccard on TSV with explicit index column name and numeric binarization threshold
python jaccard_similarity.py matrix.tsv --axis columns --index-col ProteinA --threshold 0 --out jaccard_cols.csv

# Use sparse implementation for large matrices
python jaccard_similarity.py matrix.csv --sparse --out jaccard_rows.csv

# Generate histograms for VTR-VTR, HTF-HTF, and VTR-HTF using name lists
python jaccard_similarity.py ../data/protein_interaction_matrix.csv \
  --index-col ProteinA \
  --out jaccard_rows.csv \
  --vtr-list ../data/vtrs.txt \
  --htf-list ../data/htfs.txt \
  --hist-outdir ../results/ \
  --bins 20

python jaccard_similarity.py ../data/protein_interaction_matrix.csv \
  --index-col ProteinA \
  --out jaccard_rows.csv \
  --vtr-list ../data/vtrs.txt \
  --htf-list ../data/htfs.txt \
  --hist-outdir ../results/ \
  --bins 20 --min-features 3
"""

from __future__ import annotations
import argparse
import os
import sys
from typing import Optional
import pandas as pd
import numpy as np
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import scipy.sparse as sp  # optional, only needed if --sparse
except Exception:
    sp = None


# Optional: scipy.stats for significance testing
try:
    from scipy import stats as spstats  # optional, for p-values
except Exception:
    spstats = None

# --- Mann–Whitney U test with rank-biserial effect size ---
def mannwhitney_with_effect(a, b, alternative='two-sided'):
    if spstats is None or len(a) == 0 or len(b) == 0:
        return None, None, None
    u_stat, p_val = spstats.mannwhitneyu(a, b, alternative=alternative)
    n1 = len(a); n2 = len(b)
    # rank-biserial correlation (common for MWU)
    rbc = (2.0 * u_stat) / (n1 * n2) - 1.0
    return u_stat, p_val, rbc


def autodetect_sep(path: str) -> str:
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        head = f.readline()
    # Heuristic: prefer tab if present, else comma
    if '\t' in head and head.count('\t') >= head.count(','):
        return '\t'
    return ','


def read_table(path: str, index_col: Optional[str], sep: Optional[str]) -> pd.DataFrame:
    if sep is None:
        sep = autodetect_sep(path)
    # Try reading without index first to allow user to specify index by name
    df = pd.read_csv(path, sep=sep)
    if index_col is not None:
        if index_col in df.columns:
            df = df.set_index(index_col)
        else:
            raise SystemExit(f"[error] --index-col '{index_col}' not found in columns: {list(df.columns)[:10]} ...")
    else:
        # If the first column seems non-numeric and others numeric, use it as index
        first_col = df.columns[0]
        # If first column is non-numeric OR its name suggests ID
        if not np.issubdtype(df[first_col].dtype, np.number):
            df = df.set_index(first_col)
        else:
            # If many NaNs after casting to float, assume first column was an index accidentally parsed as float
            non_numeric_count = df[first_col].isna().sum()
            if non_numeric_count > 0:
                df = pd.read_csv(path, sep=sep, index_col=0)
    return df


def binarize(df: pd.DataFrame, threshold: float) -> pd.DataFrame:
    # Convert all values to numeric (coerce errors to NaN), then fill NaN with 0
    df_num = df.apply(pd.to_numeric, errors='coerce').fillna(0.0)
    return (df_num > float(threshold))


def jaccard_dense_bool(X: np.ndarray) -> np.ndarray:
    """Row-wise Jaccard similarity for a dense boolean matrix X (n_items x n_features)."""
    X = X.astype(np.uint8, copy=False)
    # intersection via dot product in uint8 (counts of co-ones)
    inter = X @ X.T  # (n x n) counts
    row_sums = X.sum(axis=1, dtype=np.int64)  # (n,)
    # union counts = sum_i + sum_j - inter
    union = row_sums[:, None] + row_sums[None, :] - inter
    with np.errstate(divide='ignore', invalid='ignore'):
        J = inter / union
        J[union == 0] = 0.0
    return J.astype(np.float64, copy=False)


def jaccard_sparse_bool(X: sp.csr_matrix) -> np.ndarray:
    """Row-wise Jaccard similarity for a CSR boolean matrix using sparse ops."""
    if not sp.isspmatrix_csr(X):
        X = X.tocsr()
    X.data = np.ones_like(X.data)  # ensure boolean 1s
    inter = (X @ X.T).astype(np.int64)  # sparse intersections
    row_sums = np.array(X.sum(axis=1)).ravel().astype(np.int64)
    # Build union = sum_i + sum_j - inter for all pairs efficiently
    # We'll materialize dense at the end
    inter_dense = inter.toarray()
    union = row_sums[:, None] + row_sums[None, :] - inter_dense
    with np.errstate(divide='ignore', invalid='ignore'):
        J = inter_dense / union
        J[union == 0] = 0.0
    return J


def compute_jaccard(df_bool: pd.DataFrame, axis: str, use_sparse: bool) -> pd.DataFrame:
    if axis == 'columns':
        # compute column-wise Jaccard by transposing
        df_bool = df_bool.T
    labels = df_bool.index.to_list()
    if use_sparse:
        if sp is None:
            raise SystemExit("[error] scipy is required for --sparse but is not installed.")
        X = sp.csr_matrix(df_bool.values.astype(np.uint8))
        J = jaccard_sparse_bool(X)
    else:
        X = df_bool.values.astype(bool, copy=False)
        J = jaccard_dense_bool(X)
    out = pd.DataFrame(J, index=labels, columns=labels)
    return out


def read_name_list(path: str) -> list:
    """Read names from a text file. Accepts comma or whitespace separated; ignores blanks and lines starting with '#'."""
    names = []
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            for tok in re.split(r"[\s,]+", line):
                tok = tok.strip()
                if tok:
                    names.append(tok)
    # de-duplicate but preserve order
    seen = set()
    out = []
    for n in names:
        if n not in seen:
            seen.add(n)
            out.append(n)
    return out


def extract_pairs(sim: pd.DataFrame, items_a: list, items_b: list, symmetric: bool = False) -> tuple[pd.DataFrame, np.ndarray]:
    """Extract pairwise similarities between items_a and items_b from sim.
    If symmetric is True, only take the upper triangle (i<j) within items_a==items_b.
    Returns (pairs_df, values_array).
    """
    idx = sim.index
    set_idx = set(idx)
    a = [x for x in items_a if x in set_idx]
    b = [x for x in items_b if x in set_idx]
    if not a or not b:
        return pd.DataFrame(columns=['item1', 'item2', 'similarity']), np.array([], dtype=float)

    if symmetric:
        sub = sim.loc[a, a].to_numpy()
        names = a
        # upper triangle without diagonal
        triu = np.triu_indices(len(names), k=1)
        vals = sub[triu]
        pairs = [{'item1': names[i], 'item2': names[j], 'similarity': float(sub[i, j])} for i, j in zip(triu[0], triu[1])]
        return pd.DataFrame(pairs), vals
    else:
        sub = sim.loc[a, b]
        # all cross pairs
        pairs = []
        for i_name, row in sub.iterrows():
            for j_name, val in row.items():
                pairs.append({'item1': i_name, 'item2': j_name, 'similarity': float(val)})
        vals = sub.to_numpy().ravel()
        return pd.DataFrame(pairs), vals

# --- Helper to annotate pairs with feature lists and shared features ---
def add_feature_lists_and_shared(pairs_df: pd.DataFrame, df_source: pd.DataFrame) -> pd.DataFrame:
    pairs_df = pairs_df.copy()
    item1_feats = []
    item2_feats = []
    shared_col = []
    for _, row in pairs_df.iterrows():
        a = row['item1']
        b = row['item2']
        if a in df_source.index:
            a_mask = (df_source.loc[a] > 0)
            a_list = df_source.columns[a_mask].tolist()
        else:
            a_list = []
        if b in df_source.index:
            b_mask = (df_source.loc[b] > 0)
            b_list = df_source.columns[b_mask].tolist()
        else:
            b_list = []
        shared = sorted(set(a_list).intersection(b_list))
        item1_feats.append(';'.join(a_list))
        item2_feats.append(';'.join(b_list))
        shared_col.append(';'.join(shared))
    pairs_df['item1_features'] = item1_feats
    pairs_df['item2_features'] = item2_feats
    pairs_df['shared_features'] = shared_col
    return pairs_df


def plot_hist(values: np.ndarray, title: str, out_path: str, bins: int = 50) -> None:
    """Plot and save a histogram image for the given values."""
    if values.size == 0:
        return
    plt.figure()
    # Use fixed bin width of 0.1 over [0, 1]
    bin_edges = np.arange(0.0, 1.0 + 1e-9, 0.1)
    plt.hist(values, bins=bin_edges)
    plt.xlabel('Jaccard similarity')
    plt.ylabel('Count')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def write_summary(values: dict, out_path: str) -> None:
    """Write summary + tail stats for each distribution to a text file.
    values: {name: np.ndarray}
    """
    lines = []
    for name, v in values.items():
        v = v[np.isfinite(v)]
        if v.size == 0:
            lines.append(f"[{name}] no data\n")
            continue
        stats = {
            'count': v.size,
            'min': float(np.min(v)),
            'p50': float(np.percentile(v, 50)),
            'p90': float(np.percentile(v, 90)),
            'p95': float(np.percentile(v, 95)),
            'p99': float(np.percentile(v, 99)),
            'max': float(np.max(v)),
            'mean': float(np.mean(v)),
            'std': float(np.std(v, ddof=1)) if v.size > 1 else 0.0,
        }
        lines.append(f"[{name}]\n")
        for k in ['count','min','p50','p90','p95','p99','max','mean','std']:
            lines.append(f"  {k:>4}: {stats[k]}\n")
        lines.append("\n")
    with open(out_path, 'w') as f:
        f.writelines(lines)


def plot_overlaid_hist(values: dict, out_path: str, bins: int = 50) -> None:
    """Plot overlaid outline histograms as probability densities over [0,1],
    with a logarithmic y-scale to highlight long tails. The area under each
    curve integrates to ~1 (up to numerical precision).
    """
    any_nonempty = any((v.size > 0) for v in values.values())
    if not any_nonempty:
        return

    # Shared bin edges to ensure curves are comparable (fixed width 0.1)
    edges = np.arange(0.0, 1.0 + 1e-9, 0.1)

    plt.figure()

    # Small epsilon to avoid zeros on log scale (which are not plottable)
    eps = 1e-12

    for name, v in values.items():
        if v.size == 0:
            continue
        v = np.asarray(v)
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        # Manual histogram to control zero handling
        counts, _ = np.histogram(v, bins=edges, density=True)
        counts = np.maximum(counts, eps)
        centers = 0.5 * (edges[:-1] + edges[1:])
        plt.step(centers, counts, where='mid', linewidth=2, label=name)

    plt.xlabel('Jaccard similarity')
    plt.ylabel('Probability density')
    plt.title('Jaccard similarity distributions (overlaid, PDF, symlog-y)')
    plt.xlim(0.0, 1.0)

    # Use a symmetric log scale so we can include 0 as a tick (no negatives expected)
    plt.yscale('symlog', linthresh=1e-5)

    # Custom y-ticks: 0, 10^-5, 10^-4, …, 10^2
    yticks = [0.0] + [10.0 ** e for e in range(-5, 3)]  # -5, -4, ..., 2
    yticklabels = ['0'] + [rf'$10^{{{e}}}$' for e in range(-5, 3)]
    plt.yticks(yticks, yticklabels)
    plt.ylim(bottom=0)

    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def plot_cdf(values: dict, out_path: str) -> None:
    """Plot overlaid empirical CDFs for multiple groups."""
    any_nonempty = any((v.size > 0) for v in values.values())
    if not any_nonempty:
        return
    plt.figure()
    for name, v in values.items():
        v = np.asarray(v)
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        x = np.sort(v)
        y = np.arange(1, x.size + 1) / x.size
        plt.plot(x, y, label=name, linewidth=2)
    plt.xlabel('Jaccard similarity')
    plt.ylabel('CDF (fraction ≤ x)')
    plt.title('Jaccard similarity CDFs (overlaid)')
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

# --- Violin plot helper ---
def plot_violin(values: dict, out_path: str) -> None:
    """Plot a violin plot for multiple groups over [0,1]."""
    # Collect non-empty arrays and labels in a stable order
    data = []
    labels = []
    for name, v in values.items():
        v = np.asarray(v)
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        data.append(v)
        labels.append(name)
    if not data:
        return

    plt.figure()
    parts = plt.violinplot(data, showmeans=True, showmedians=True, showextrema=True)
    # Slight transparency so overlaps/quantiles are visible
    for pc in parts.get('bodies', []):
        pc.set_alpha(0.6)
    plt.xticks(np.arange(1, len(labels) + 1), labels, rotation=0)
    plt.ylim(0.0, 1.0)
    plt.ylabel('Jaccard similarity')
    plt.title('Jaccard similarity distributions (violin)')
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


# --- CCDF plotting helper ---
def plot_ccdf(values: dict, out_path: str) -> None:
    """Plot overlaid empirical CCDFs (1 - CDF) with log-y to emphasize tails."""
    any_nonempty = any((v.size > 0) for v in values.values())
    if not any_nonempty:
        return
    plt.figure()
    for name, v in values.items():
        v = np.asarray(v)
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        x = np.sort(v)
        y = 1.0 - (np.arange(1, x.size + 1) / x.size)
        plt.plot(x, y, label=name, linewidth=2)
    plt.xlabel('Jaccard similarity')
    plt.ylabel('1 - CDF (fraction ≥ x)')
    plt.title('Jaccard similarity CCDFs (overlaid, log-y)')
    plt.xlim(0.0, 1.0)
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()




def main():
    p = argparse.ArgumentParser(description="Compute Jaccard similarity matrix from a feature table.")
    p.add_argument('table', help='Input table (CSV/TSV). Rows = items; columns = features.')
    p.add_argument('--index-col', help='Name of the column to use as row index (ID). If omitted, auto-detect.', default=None)
    p.add_argument('--sep', help='Field delimiter: "," or "\t". Default: auto-detect.', default=None)
    p.add_argument('--threshold', type=float, default=0.0, help='Binarization threshold: value > threshold => 1. Default: 0.')
    p.add_argument('--axis', choices=['rows', 'columns'], default='rows', help='Compute Jaccard across rows or columns. Default: rows')
    p.add_argument('--sparse', action='store_true', help='Use sparse implementation (recommended for very wide/tall matrices).')
    p.add_argument('--out', required=True, help='Output CSV path for the Jaccard similarity matrix.')
    p.add_argument('--vtr-list', help='Path to text file with VTR names (one per line, or comma/whitespace separated).', default=None)
    p.add_argument('--htf-list', help='Path to text file with HTF names (one per line, or comma/whitespace separated).', default=None)
    p.add_argument('--hist-outdir', help='If provided with --vtr-list and --htf-list, write histograms and per-pair CSVs here.', default=None)
    p.add_argument('--bins', type=int, default=50, help='Number of bins for histograms (default: 50).')
    p.add_argument('--ccdf', action='store_true', help='Also write 1-CDF (CCDF) plot with log-y to emphasize the upper tail.')
    p.add_argument('--tf-cols', default='../data/tf.txt',
                   help='Path to a list of TF column names to subset for TF-only Jaccard (default: ../data/tf.txt).')
    p.add_argument('--cof-cols', default='../data/cof.txt',
                   help='Path to a list of cofactor column names to subset for Cof-only Jaccard (default: ../data/cof.txt).')
    p.add_argument('--min-features', type=int, default=0,
                   help='If >0, drop items (rows) that have fewer than this many 1s after binarization. Applied before Jaccard.')
    args = p.parse_args()

    df = read_table(args.table, index_col=args.index_col, sep=args.sep)
    if df.empty:
        raise SystemExit('[error] Input table contains no data after parsing.')

    # Binarize values
    df_bool = binarize(df, threshold=args.threshold)

    # Filter rows with too few features if requested
    if args.min_features and args.min_features > 0:
        row_counts = df_bool.sum(axis=1)
        keep_mask = row_counts >= args.min_features
        dropped = (~keep_mask).sum()
        if dropped > 0:
            print(f"[info] Dropping {dropped} items with < {args.min_features} features", file=sys.stderr)
        df_bool = df_bool.loc[keep_mask]

    # Compute Jaccard
    sim = compute_jaccard(df_bool, axis=args.axis, use_sparse=args.sparse)

    # Read TF/CoF column name lists if present
    tf_list = read_name_list(args.tf_cols) if os.path.exists(args.tf_cols) else []
    cof_list = read_name_list(args.cof_cols) if os.path.exists(args.cof_cols) else []

    # Write output
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    sim.to_csv(args.out)
    print(f"[ok] Wrote Jaccard {args.axis} similarity matrix: {args.out}")

    # Optionally generate histograms for VTR-VTR, HTF-HTF, and VTR-HTF pairs
    if args.vtr_list and args.htf_list and args.hist_outdir:
        os.makedirs(args.hist_outdir, exist_ok=True)

        vtrs = read_name_list(args.vtr_list)
        htfs = read_name_list(args.htf_list)

        # Build VTR virus map
        def _vtr_virus_map(vtr_names):
            virus = {}
            for name in vtr_names:
                if '_' in name:
                    virus[name] = name.split('_', 1)[1]
                else:
                    virus[name] = 'UNK'
            return virus
        vtr_virus = _vtr_virus_map(vtrs)

        # Report any missing names
        idx_set = set(sim.index)
        missing_v = sorted(set(vtrs) - idx_set)
        missing_h = sorted(set(htfs) - idx_set)
        if missing_v:
            print(f"[warn] {len(missing_v)} VTR names not found in similarity matrix index (showing up to 10): {missing_v[:10]}", file=sys.stderr)
        if missing_h:
            print(f"[warn] {len(missing_h)} HTF names not found in similarity matrix index (showing up to 10): {missing_h[:10]}", file=sys.stderr)

        # VTR-VTR
        df_vv, vals_vv = extract_pairs(sim, vtrs, vtrs, symmetric=True)
        df_vv = add_feature_lists_and_shared(df_vv, df_bool)
        vv_csv = os.path.join(args.hist_outdir, 'pairs_vtr_vtr.csv')
        df_vv.to_csv(vv_csv, index=False)
        plot_hist(vals_vv, 'Jaccard similarity: VTR vs VTR', os.path.join(args.hist_outdir, 'hist_vtr_vtr.png'), bins=args.bins)

        # HTF-HTF
        df_hh, vals_hh = extract_pairs(sim, htfs, htfs, symmetric=True)
        df_hh = add_feature_lists_and_shared(df_hh, df_bool)
        hh_csv = os.path.join(args.hist_outdir, 'pairs_htf_htf.csv')
        df_hh.to_csv(hh_csv, index=False)
        plot_hist(vals_hh, 'Jaccard similarity: HTF vs HTF', os.path.join(args.hist_outdir, 'hist_htf_htf.png'), bins=args.bins)

        # VTR-HTF
        df_vh, vals_vh = extract_pairs(sim, vtrs, htfs, symmetric=False)
        df_vh = add_feature_lists_and_shared(df_vh, df_bool)
        vh_csv = os.path.join(args.hist_outdir, 'pairs_vtr_htf.csv')
        df_vh.to_csv(vh_csv, index=False)
        plot_hist(vals_vh, 'Jaccard similarity: VTR vs HTF', os.path.join(args.hist_outdir, 'hist_vtr_htf.png'), bins=args.bins)

        # Overlaid histograms (all three together)
        grouped = {
            'VTR–VTR': vals_vv,
            'HTF–HTF': vals_hh,
            'VTR–HTF': vals_vh,
        }
        plot_overlaid_hist(grouped, os.path.join(args.hist_outdir, 'hist_overlaid.pdf'), bins=args.bins)

        # Also plot same-virus vs different-virus for ALL columns (VTR–VTR pairs)
        same_vals_all = []
        diff_vals_all = []
        for _, r in df_vv.iterrows():
            a = r['item1']; b = r['item2']
            va = vtr_virus.get(a, 'UNK'); vb = vtr_virus.get(b, 'UNK')
            if va == vb:
                same_vals_all.append(r['similarity'])
            else:
                diff_vals_all.append(r['similarity'])
        grouped_vtr_byvirus_all = {
            'VTR–VTR (same virus)': np.array(same_vals_all, dtype=float),
            'VTR–VTR (different virus)': np.array(diff_vals_all, dtype=float),
        }
        plot_overlaid_hist(grouped_vtr_byvirus_all, os.path.join(args.hist_outdir, 'hist_overlaid_vtr_byvirus.pdf'), bins=args.bins)

        # Compute significance for ALL columns (same vs different virus)
        u_stat, p_val, rbc = mannwhitney_with_effect(same_vals_all, diff_vals_all)
        with open(os.path.join(args.hist_outdir, 'pvals_vtr_byvirus_all.txt'), 'w') as f:
            if u_stat is not None:
                f.write(f"Mann-Whitney U (same vs diff virus, ALL): U={u_stat}, p={p_val}, rbc={rbc}\n")
            else:
                f.write("scipy.stats not available or one of the groups is empty; p-value not computed.\n")

        # VTR–VTR vs VTR–HTF distributions (ALL columns)
        u_stat2, p_val2, rbc2 = mannwhitney_with_effect(vals_vv, vals_vh)
        with open(os.path.join(args.hist_outdir, 'pvals_vtr_vs_htf_all.txt'), 'w') as f:
            if u_stat2 is not None:
                f.write(f"Mann-Whitney U (VTR–VTR vs VTR–HTF, ALL): U={u_stat2}, p={p_val2}, rbc={rbc2}\n")
            else:
                f.write("scipy.stats not available or one of the groups is empty; p-value not computed.\n")

        # Helper to write same/diff virus VTR–VTR pair CSVs for a tag
        def _write_subset_vtr_virus_tables(tag, df_vv_s, vtr_virus_map, df_source):
            df_vv_s = add_feature_lists_and_shared(df_vv_s, df_source)
            same_rows = []
            diff_rows = []
            for _, row in df_vv_s.iterrows():
                a = row['item1']; b = row['item2']
                va = vtr_virus_map.get(a, 'UNK'); vb = vtr_virus_map.get(b, 'UNK')
                if va == vb:
                    same_rows.append(row)
                else:
                    diff_rows.append(row)
            same_df = pd.DataFrame(same_rows)
            diff_df = pd.DataFrame(diff_rows)
            same_csv = os.path.join(args.hist_outdir, f'pairs_vtr_vtr_samevirus_{tag}.csv')
            diff_csv = os.path.join(args.hist_outdir, f'pairs_vtr_vtr_diffvirus_{tag}.csv')
            same_df.to_csv(same_csv, index=False)
            diff_df.to_csv(diff_csv, index=False)

        # --- Additional overlaid histograms by column subsets: TF, Cofactor, Other ---
        def _subset_and_plot(cols_subset, tag):
            # Intersect with available columns
            use_cols = [c for c in df_bool.columns if c in set(cols_subset)]
            if len(use_cols) == 0:
                print(f"[warn] No columns found for subset '{tag}'; skipping.", file=sys.stderr)
                return
            df_sub = df_bool.loc[:, use_cols]
            sim_sub = compute_jaccard(df_sub, axis=args.axis, use_sparse=args.sparse)
            # VTR–VTR pairs
            df_vv_s, vals_vv_s = extract_pairs(sim_sub, vtrs, vtrs, symmetric=True)
            # HTF–HTF pairs
            df_hh_s, vals_hh_s = extract_pairs(sim_sub, htfs, htfs, symmetric=True)
            # VTR–HTF pairs and write CSV (move this up for later use)
            df_vh_s, _ = extract_pairs(sim_sub, vtrs, htfs, symmetric=False)
            df_vh_s = add_feature_lists_and_shared(df_vh_s, df_sub)
            vh_csv_s = os.path.join(args.hist_outdir, f'pairs_vtr_htf_{tag}.csv')
            df_vh_s.to_csv(vh_csv_s, index=False)
            # Split VTR–VTR pairs by virus
            same_vals = []
            diff_vals = []
            for _, r in df_vv_s.iterrows():
                a = r['item1']; b = r['item2']
                va = vtr_virus.get(a, 'UNK'); vb = vtr_virus.get(b, 'UNK')
                if va == vb:
                    same_vals.append(r['similarity'])
                else:
                    diff_vals.append(r['similarity'])
            grouped_s = {
                'VTR–VTR (same virus)': np.array(same_vals, dtype=float),
                'VTR–VTR (different virus)': np.array(diff_vals, dtype=float),
            }
            _write_subset_vtr_virus_tables(tag, df_vv_s, vtr_virus, df_sub)
            plot_overlaid_hist(grouped_s, os.path.join(args.hist_outdir, f'hist_overlaid_{tag}.pdf'), bins=args.bins)
            # Significance test for this subset (same vs different virus)
            pval_path = os.path.join(args.hist_outdir, f'pvals_vtr_byvirus_{tag}.txt')
            u_stat_s, p_val_s, rbc_s = mannwhitney_with_effect(same_vals, diff_vals)
            with open(pval_path, 'w') as f:
                if u_stat_s is not None:
                    f.write(f"Mann-Whitney U (same vs diff virus, {tag}): U={u_stat_s}, p={p_val_s}, rbc={rbc_s}\n")
                else:
                    f.write("scipy.stats not available or one of the groups is empty; p-value not computed.\n")
            # Now plot 3-way overlaid hist for this subset: VTR–VTR, HTF–HTF, VTR–HTF (mirroring original grouped)
            grouped_all = {
                'VTR–VTR': vals_vv_s,
                'HTF–HTF': vals_hh_s,
                'VTR–HTF': df_vh_s['similarity'].to_numpy() if not df_vh_s.empty else np.array([], dtype=float),
            }
            plot_overlaid_hist(grouped_all, os.path.join(args.hist_outdir, f'hist_overlaid_{tag}_vtr_htf.pdf'), bins=args.bins)

            # Mann–Whitney U test for VTR–VTR vs VTR–HTF for this subset
            pval_path2 = os.path.join(args.hist_outdir, f'pvals_vtr_vs_htf_{tag}.txt')
            u_stat_vh, p_val_vh, rbc_vh = mannwhitney_with_effect(
                vals_vv_s,
                df_vh_s['similarity'].to_numpy() if not df_vh_s.empty else np.array([], dtype=float)
            )
            with open(pval_path2, 'w') as f:
                if u_stat_vh is not None:
                    f.write(f"Mann-Whitney U (VTR–VTR vs VTR–HTF, {tag}): U={u_stat_vh}, p={p_val_vh}, rbc={rbc_vh}\n")
                else:
                    f.write("scipy.stats not available or one of the groups is empty; p-value not computed.\n")

        # TF-only plot
        if tf_list:
            _subset_and_plot(tf_list, 'tf')
        else:
            print(f"[warn] TF list not found or empty: {args.tf_cols}", file=sys.stderr)

        # Cofactor-only plot
        if cof_list:
            _subset_and_plot(cof_list, 'cof')
        else:
            print(f"[warn] Cof list not found or empty: {args.cof_cols}", file=sys.stderr)

        # Other-columns plot (not in TF or Cof lists)
        if tf_list or cof_list:
            tfcof = set(tf_list) | set(cof_list)
            other_cols = [c for c in df_bool.columns if c not in tfcof]
            if other_cols:
                _subset_and_plot(other_cols, 'other')
            else:
                print("[warn] No 'other' columns (everything is in TF or Cof).", file=sys.stderr)

        # Overlaid CDFs
        plot_cdf(grouped, os.path.join(args.hist_outdir, 'cdf_overlaid.png'))
        plot_violin(grouped, os.path.join(args.hist_outdir, 'violin_overlaid.png'))
        if args.ccdf:
            plot_ccdf(grouped, os.path.join(args.hist_outdir, 'ccdf_overlaid.png'))

        # Write numeric summary (includes max/min and percentiles up to 99th)
        write_summary(grouped, os.path.join(args.hist_outdir, 'summary_stats.txt'))

        print(f"[ok] Wrote pair tables to: {args.hist_outdir}\n      - VTR-VTR pairs: {len(df_vv)}\n      - HTF-HTF pairs: {len(df_hh)}\n      - VTR-HTF pairs: {len(df_vh)}\n[ok] Wrote histograms: hist_vtr_vtr.png, hist_htf_htf.png, hist_vtr_htf.png")
    elif args.hist_outdir and not (args.vtr_list and args.htf_list):
        print('[warn] --hist-outdir was provided but missing --vtr-list/--htf-list; skipping histogram generation.', file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)