#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
viz_wgcna_se2_fast.py
ויזואליזציה מהירה ואינפורמטיבית לקבצי WGCNA↔SpeakEasy2 הקיימים.

ייבוא קבצים (בשם קבוע) מתוך --indir, ויצירת תמונות PNG לתיקיית --outdir:
  - heatmap_jaccard_top.png       : Heatmap של Jaccard (Top rows/cols), עם נקודות לסיג' (q<qthr)
  - heatmap_significant_top.png   : Heatmap בינארי של תאים מובהקים (q<qthr) (Top rows/cols)
  - bar_top_pairs_jaccard.png     : טופ 25 התאמות Hungarian לפי Jaccard
  - hist_purity.png               : היסט' Purity של מודולים (לבן הזוג שנבחר)
  - hist_f1.png                   : היסט' F1 של מודולים
  - hist_max_overlap_se2.png      : התפלגות מקס' Overlap לכל קהילת SE2 מול מודולים
  - hist_max_overlap_wgcna.png    : התפלגות מקס' Overlap לכל מודול WGCNA מול קהילות
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read_matrix_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    # אם העמודה הראשונה נראית כמו שמות שורות, קבע כאינדקס
    if df.shape[1] > 1 and not pd.api.types.is_numeric_dtype(df.iloc[:,0]):
        df = df.set_index(df.columns[0])
    # המרה לנומרי
    df = df.apply(pd.to_numeric, errors="coerce")
    return df

def read_scores_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df

def select_top(J: pd.DataFrame, top_rows: int, top_cols: int) -> pd.DataFrame:
    # בחירת שורות/עמודות עם המקסימום הגבוה ביותר (לפי Jaccard)
    row_order = J.max(axis=1).sort_values(ascending=False).head(top_rows).index
    col_order = J.max(axis=0).sort_values(ascending=False).head(top_cols).index
    return J.loc[row_order, col_order]

def plot_heatmap(mat: pd.DataFrame, out_png: str, title: str):
    fig, ax = plt.subplots(figsize=(max(6, 0.25*mat.shape[1]), max(6, 0.25*mat.shape[0])))
    im = ax.imshow(mat.values, aspect="auto")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=90)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("value")
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

def plot_heatmap_with_sig(J: pd.DataFrame, Q: pd.DataFrame, out_png: str, qthr: float, title: str):
    # נחתוך את Q לאותו סדר וגודל
    Qs = Q.reindex(index=J.index, columns=J.columns)
    sig = (Qs < qthr)
    fig, ax = plt.subplots(figsize=(max(6, 0.25*J.shape[1]), max(6, 0.25*J.shape[0])))
    im = ax.imshow(J.values, aspect="auto")
    # נקודות איפה שיש מובהקות
    ys, xs = np.where(sig.values)
    ax.scatter(xs, ys, s=6)  # לא להגדיר צבעים — ברירת מחדל
    ax.set_xticks(range(J.shape[1])); ax.set_xticklabels(J.columns, rotation=90)
    ax.set_yticks(range(J.shape[0])); ax.set_yticklabels(J.index)
    ax.set_title(title + f" (dots: q<{qthr})")
    cbar = fig.colorbar(im, ax=ax); cbar.set_label("Jaccard")
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

def plot_binary_heatmap(B: pd.DataFrame, out_png: str, title: str):
    fig, ax = plt.subplots(figsize=(max(6, 0.25*B.shape[1]), max(6, 0.25*B.shape[0])))
    im = ax.imshow(B.values.astype(float), aspect="auto")
    ax.set_xticks(range(B.shape[1])); ax.set_xticklabels(B.columns, rotation=90)
    ax.set_yticks(range(B.shape[0])); ax.set_yticklabels(B.index)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax); cbar.set_label("binary")
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

def plot_hist(series: pd.Series, out_png: str, title: str, xlabel: str):
    s = pd.to_numeric(series, errors="coerce").dropna()
    fig, ax = plt.subplots(figsize=(7,5))
    ax.hist(s.values, bins=30)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("count")
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

def plot_bar(labels: pd.Series, values: pd.Series, out_png: str, title: str, xlabel: str):
    # תוויות ארוכות: נגביל ל־25 ברירות־מחדל
    n = min(25, len(values))
    labs = labels.iloc[:n]
    vals = values.iloc[:n]
    fig, ax = plt.subplots(figsize=(max(8, 0.4*n), 6))
    ax.bar(range(n), vals.values)
    ax.set_xticks(range(n))
    ax.set_xticklabels(labs.values, rotation=90)
    ax.set_title(title)
    ax.set_ylabel(xlabel)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, help="Folder containing the CSVs")
    ap.add_argument("--outdir", required=True, help="Where to write PNGs")
    ap.add_argument("--top-rows", type=int, default=40)
    ap.add_argument("--top-cols", type=int, default=40)
    ap.add_argument("--qthr", type=float, default=0.05)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    p = lambda f: os.path.join(args.indir, f)

    # --- Read matrices ---
    J = read_matrix_csv(p("jaccard.csv"))
    O = read_matrix_csv(p("overlap.csv"))
    Q = read_matrix_csv(p("qvalues_fdr.csv"))

    # היישור: שמות השורות/עמודות צריכים להתאים ביניהם; אם יש אי־התאמה, נחיתוך ל־intersection
    common_rows = J.index.intersection(Q.index)
    common_cols = J.columns.intersection(Q.columns)
    J = J.loc[common_rows, common_cols]
    Q = Q.loc[common_rows, common_cols]

    # בחירת טופ־שורות/עמודות
    J_top = select_top(J, args.top_rows, args.top_cols)
    Q_top = Q.loc[J_top.index, J_top.columns]

    # Heatmaps
    plot_heatmap_with_sig(J_top, Q_top,
                          out_png=os.path.join(args.outdir, "heatmap_jaccard_top.png"),
                          qthr=args.qthr,
                          title="Top Jaccard (rows/cols by max J)")

    # Heatmap בינארי של מובהקות
    B = (Q_top < args.qthr).astype(int)
    plot_binary_heatmap(B,
                        out_png=os.path.join(args.outdir, "heatmap_significant_top.png"),
                        title=f"Significant cells (q<{args.qthr}) - top rows/cols")

    # --- Best pairs bar chart ---
    bp = pd.read_csv(p("best_pairs_hungarian.csv"))
    bp = bp.sort_values("jaccard", ascending=False).reset_index(drop=True)
    labels = bp.apply(lambda r: f"WGCNA:{r['wgcna']} → SE2:{r['se2']}", axis=1)
    plot_bar(labels, bp["jaccard"],
             out_png=os.path.join(args.outdir, "bar_top_pairs_jaccard.png"),
             title="Top Hungarian matches by Jaccard",
             xlabel="Jaccard")

    # --- Per-module distributions ---
    scores = pd.read_csv(p("per_module_scores.csv"))
    plot_hist(scores["purity"],
              out_png=os.path.join(args.outdir, "hist_purity.png"),
              title="Module purity distribution (WGCNA → matched SE2)",
              xlabel="purity")
    plot_hist(scores["F1"],
              out_png=os.path.join(args.outdir, "hist_f1.png"),
              title="Module F1 distribution (WGCNA ↔ matched SE2)",
              xlabel="F1")

    # --- Overlap summaries ---
    # מקסימום לכל קהילה/מודול (כמה הקבוצה הקטנה מוכללת)
    # וודא ש־O מיושר כמו J/Q (אם שמות לא זהים)
    if not O.index.equals(J.index) or not O.columns.equals(J.columns):
        rows = O.index.intersection(J.index)
        cols = O.columns.intersection(J.columns)
        O = O.loc[rows, cols]

    se2_max = O.max(axis=0)
    wgcna_max = O.max(axis=1)
    plot_hist(se2_max,
              out_png=os.path.join(args.outdir, "hist_max_overlap_se2.png"),
              title="Max Overlap per SE2 cluster vs WGCNA modules",
              xlabel="max overlap")
    plot_hist(wgcna_max,
              out_png=os.path.join(args.outdir, "hist_max_overlap_wgcna.png"),
              title="Max Overlap per WGCNA module vs SE2 clusters",
              xlabel="max overlap")

    # הודעת הצלחה קצרה
    print("Wrote figures to:", args.outdir)

if __name__ == "__main__":
    main()
