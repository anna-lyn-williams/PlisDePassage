#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore
from matplotlib.colors import Normalize


# ----------------------------
# Label cleanup
# ----------------------------
def clean_labels(labels):
    return [str(x).replace("Area", "") for x in labels]


# ----------------------------
# Ordering (optional but preserved)
# ----------------------------
DESIRED_CYTO_ORDER = [
    "AreahIP7.IPS",
    "AreahIP4.IPS",
    "AreahIP5.IPS",
    "AreahIP8.IPS",
    "AreahIP6.IPS",
    "Area7P.SPL",
    "Area7A.SPL",
]

DESIRED_ROI_ORDER = [
    "V3ab_IPS0",
    "IPS0_IPS1",
    "IPS1_IPS2",
    "IPS2_IPS3",
]

# Heatmap colorbar limits
CBAR_MIN = 0.0
CBAR_MAX = 1.0

def apply_cyto_order(mat):
    present = [c for c in DESIRED_CYTO_ORDER if c in mat.index]
    extras = [c for c in mat.index if c not in DESIRED_CYTO_ORDER]
    return mat.reindex(present + extras)


def apply_roi_order(mat):
    present = [c for c in DESIRED_ROI_ORDER if c in mat.columns]
    extras = [c for c in mat.columns if c not in DESIRED_ROI_ORDER]
    return mat.reindex(columns=present + extras)


# ----------------------------
# Data processing
# ----------------------------
def compute_overlap(df):
    df = df.copy()
    df["hit"] = pd.to_numeric(df["hit"], errors="coerce").fillna(0).astype(int)

    roi_present = (
        df.groupby(["subject", "hemi", "roi"])["vertex"]
        .nunique()
        .reset_index(name="n_vertices")
    )
    roi_present["roi_present"] = roi_present["n_vertices"] > 0

    subject_hits = (
        df.groupby(["subject", "hemi", "roi", "hit_labels"])["hit"]
        .max()
        .reset_index()
    )

    subject_hits = subject_hits.merge(
        roi_present[["subject", "hemi", "roi", "roi_present"]],
        on=["subject", "hemi", "roi"],
        how="left",
    )
    subject_hits = subject_hits[subject_hits["roi_present"].fillna(False)]

    # Gate the denominator on cyto label presence too, so this matches the D_
    # pipeline: % is computed only over subjects who actually HAVE that cyto
    # label in that hemisphere. Subjects with the label absent are excluded
    # (not counted as non-overlap, which would deflate the percentage).
    if "label_present" in df.columns:
        df["label_present"] = (
            pd.to_numeric(df["label_present"], errors="coerce").fillna(0).astype(int)
        )
        cyto_present = (
            df.groupby(["subject", "hemi", "hit_labels"])["label_present"]
            .max()
            .reset_index(name="cyto_present")
        )
        cyto_present["cyto_present"] = cyto_present["cyto_present"] > 0

        subject_hits = subject_hits.merge(
            cyto_present,
            on=["subject", "hemi", "hit_labels"],
            how="left",
        )
        subject_hits["cyto_present"] = subject_hits["cyto_present"].fillna(False)
        subject_hits = subject_hits[subject_hits["cyto_present"]]

    return (
        subject_hits.groupby(["hemi", "hit_labels", "roi"])["hit"]
        .mean()
        .reset_index(name="subject_percentage")
    )


def make_matrix(overlap, hemi):
    mat = overlap[overlap["hemi"] == hemi].pivot(
        index="hit_labels", columns="roi", values="subject_percentage"
    )
    mat = apply_cyto_order(mat)
    mat = apply_roi_order(mat)

    # clean labels
    mat.index = clean_labels(mat.index)
    mat.columns = clean_labels(mat.columns)

    return mat


# ----------------------------
# Clustering
# ----------------------------
def cluster(mat, metric="correlation", method="average", min_shared=2):
    """
    Hierarchical clustering of cyto maps (rows).

    Returns (row-reordered matrix, linkage Z, orig_labels).

    For metric="correlation", uses PAIRWISE DELETION: each pair's correlation
    is computed only over ROIs where both rows have real (non-NaN) data. No
    values are imputed. If fewer than `min_shared` shared ROIs, or either row
    is constant over the shared ROIs, the distance is set to 1.0 (far).

    orig_labels correspond to the ORIGINAL row order used to build Z, and must
    be the labels passed to dendrogram() to avoid a label/leaf mismatch.
    """
    mat = mat.loc[~mat.isna().all(axis=1)]
    if mat.shape[0] < 2:
        raise ValueError("Need at least 2 cyto maps with data to cluster.")

    orig_labels = list(mat.index)

    if metric != "correlation":
        X = mat.to_numpy(dtype=float)
        Z = linkage(pdist(X, metric=metric), method=method)
        leaves = dendrogram(Z, no_plot=True)["leaves"]
        return mat.iloc[leaves], Z, orig_labels

    X = mat.to_numpy(dtype=float)
    n = X.shape[0]
    D = np.zeros((n, n), dtype=float)

    for i in range(n):
        xi = X[i]
        for j in range(i + 1, n):
            xj = X[j]
            m = ~np.isnan(xi) & ~np.isnan(xj)
            k = int(m.sum())

            if k < min_shared:
                d = 1.0
            else:
                a, b = xi[m], xj[m]
                if np.std(a) == 0 or np.std(b) == 0:
                    d = 1.0
                else:
                    r = np.corrcoef(a, b)[0, 1]
                    d = 1.0 if np.isnan(r) else float(np.clip(1.0 - r, 0.0, 2.0))

            D[i, j] = D[j, i] = d

    Z = linkage(squareform(D, checks=False), method=method)
    leaves = dendrogram(Z, no_plot=True)["leaves"]
    return mat.iloc[leaves], Z, orig_labels


# ----------------------------
# Visualization settings
# ----------------------------
plt.rcParams.update({
    "font.size": 15,
    "lines.linewidth": 3,  
})

def plot_heatmap(mat, outpath, title):
    outpath = Path(outpath).with_suffix(".pdf")
    outpath.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 8))  # square figure

    norm = Normalize(vmin=CBAR_MIN, vmax=CBAR_MAX)

    cmap = sns.color_palette("rocket", as_cmap=True)

    sns.heatmap(
        mat,
        annot=True,
        fmt=".2f",
        mask=mat.isna(),
        cmap=cmap,                
        norm=norm,                 
        cbar_kws={
            "label": "% overlap",
            "ticks": [0, 0.25, 0.5, 0.75, 1.0],
        },
        square=False,
        ax=ax,
    )

    ax.set_title(title)

    plt.savefig(
        outpath,
        format="pdf",              # vector output
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.close()


def plot_dendrogram(Z, labels, outpath, title):
    outpath = Path(outpath).with_suffix(".pdf")
    outpath.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))

    dendrogram(
        Z,
        labels=clean_labels(labels),
        orientation="right",
        leaf_font_size=15,
        color_threshold=0,
        above_threshold_color="black",
        ax=ax,
    )

    # Dendrogram lines
    for line in ax.collections:
        line.set_color("black")
        line.set_linewidth(3)

    # Axis styling: box off, thick axes
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(2)
    ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(axis="both", width=2)

    ax.set_title(title)

    plt.savefig(
        outpath,
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.close()

def make_matrix_avg(overlap):
    mat = overlap.pivot(
        index="hit_labels",
        columns="roi",
        values="subject_percentage"
    )

    mat = apply_cyto_order(mat)
    mat = apply_roi_order(mat)

    mat.index = clean_labels(mat.index)
    mat.columns = clean_labels(mat.columns)

    return mat
# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--metric", default="correlation")
    ap.add_argument("--method", default="average")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.csv)
    overlap = compute_overlap(df)

    for hemi in sorted(overlap["hemi"].unique()):
        mat = make_matrix(overlap, hemi)

        # heatmap
        plot_heatmap(
            mat,
            outdir / f"heatmap_{hemi}.png",
            f"Cyto × ROI overlap — {hemi.upper()}",
        )

        # clustering
        mat_ord, Z, orig_labels = cluster(mat, metric=args.metric, method=args.method)

        # dendrogram (use orig_labels: Z was built in original row order)
        plot_dendrogram(
            Z,
            orig_labels,
            outdir / f"dendrogram_{hemi}.png",
            f"Cyto clustering — {hemi.upper()}",
        )

        # clustered heatmap
        plot_heatmap(
            mat_ord,
            outdir / f"heatmap_clustered_{hemi}.png",
            f"Clustered Cyto × ROI — {hemi.upper()}",
        )

    # ============================
    # LH + RH AVERAGE
    # ============================
    overlap_avg = (
        overlap
        .groupby(["hit_labels", "roi"], as_index=False)["subject_percentage"]
        .mean()
    )

    mat_avg = make_matrix_avg(overlap_avg)

    # heatmap (average)
    plot_heatmap(
        mat_avg,
        outdir / "heatmap_avg_LH_RH.pdf",
        "Cyto × ROI overlap — LH + RH average",
    )

    # clustering (average)
    mat_avg_ord, Z_avg, orig_labels_avg = cluster(
        mat_avg,
        metric=args.metric,
        method=args.method,
    )

    # dendrogram (average) — use orig_labels_avg to match Z's row order
    plot_dendrogram(
        Z_avg,
        orig_labels_avg,
        outdir / "dendrogram_avg_LH_RH.pdf",
        "Cyto clustering — LH + RH average",
    )

    # clustered heatmap (average)
    plot_heatmap(
        mat_avg_ord,
        outdir / "heatmap_clustered_avg_LH_RH.pdf",
        "Clustered Cyto × ROI — LH + RH average",
    )
    print(f"Done → {outdir}")


if __name__ == "__main__":
    main()