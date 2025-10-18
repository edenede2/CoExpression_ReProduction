#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compare WGCNA modules vs SpeakEasy2 clusters (per-gene membership tables).
Outputs CSVs + JSON summary + matplotlib figures (PNG/PDF) + README.

Usage:
  python compare_wgcna_se2.py --wgcna <path> --se2 <path> --outdir <path>
  [--gene-col-wgcna COL] [--module-col COL] [--gene-col-se2 COL] [--se2-col COL]
  [--drop-grey] [--max-pairs-per-module INT] [--n-perm INT]

Typical for your files:
  --gene-col-wgcna "Gene Symbol" --module-col "Cluster ID"
  --gene-col-se2 "Gene.ID"       --se2-col "Cluster.ID"
"""

import argparse
import os
import sys
import json
from typing import List, Optional, Dict
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import matthews_corrcoef, average_precision_score, precision_recall_curve
from itertools import combinations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- IO helpers ----------
def smart_read(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep=None, engine="python")
    except Exception:
        try:
            return pd.read_csv(path, sep="\t")
        except Exception:
            return pd.read_csv(path)

def auto_find_column(columns: List[str], candidates: List[str]) -> Optional[str]:
    lower = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    for c in columns:
        cl = c.lower()
        for cand in candidates:
            if cand.lower() in cl:
                return c
    return None

def standardize(wgcna_df: pd.DataFrame, se2_df: pd.DataFrame,
                gene_col_w: Optional[str], wgcna_col: Optional[str],
                gene_col_s: Optional[str], se2_col: Optional[str],
                drop_grey: bool) -> pd.DataFrame:
    if gene_col_w is None:
        gene_col_w = auto_find_column(wgcna_df.columns.tolist(),
                                      ["Gene.ID","Gene Symbol","Gene.Symbol","gene","symbol","ensembl"])
    if wgcna_col is None:
        wgcna_col = auto_find_column(wgcna_df.columns.tolist(),
                                     ["Cluster.ID","Cluster ID","module","cluster","community","modulecolor","module_color"])
    if gene_col_s is None:
        gene_col_s = auto_find_column(se2_df.columns.tolist(),
                                      ["Gene.ID","Gene.Symbol","Gene Symbol","gene","symbol","ensembl"])
    if se2_col is None:
        se2_col = auto_find_column(se2_df.columns.tolist(),
                                   ["Cluster.ID","Cluster ID","cluster","community","label"])

    missing = []
    if gene_col_w is None: missing.append("gene_col_wgcna")
    if wgcna_col is None:  missing.append("module_col (WGCNA)")
    if gene_col_s is None: missing.append("gene_col_se2")
    if se2_col is None:    missing.append("se2_col")
    if missing:
        raise ValueError(f"Could not auto-detect: {missing}\n"
                         f"WGCNA cols: {list(wgcna_df.columns)}\nSE2 cols: {list(se2_df.columns)}")

    # Include tissue in the analysis to handle genes appearing in multiple tissues
    ww = wgcna_df[[gene_col_w, wgcna_col, 'Tissue']].copy()
    ww.columns = ["gene","wgcna","tissue"]
    ss = se2_df[[gene_col_s, se2_col, 'Tissue']].copy()
    ss.columns = ["gene","se2","tissue"]

    for c in ["gene","wgcna","tissue"]:
        ww[c] = ww[c].astype(str).str.strip()
    for c in ["gene","se2","tissue"]:
        ss[c] = ss[c].astype(str).str.strip()
    
    # Create unique gene-tissue identifiers
    ww["gene"] = ww["gene"] + "_" + ww["tissue"]
    ss["gene"] = ss["gene"] + "_" + ss["tissue"]
    
    # Drop the separate tissue column as it's now incorporated into gene
    ww = ww.drop("tissue", axis=1)
    ss = ss.drop("tissue", axis=1)

    ww = ww.dropna().query("gene != '' and wgcna != ''").copy()
    ss = ss.dropna().query("gene != '' and se2 != ''").copy()

    if drop_grey:
        ww = ww[ww["wgcna"].str.lower()!="grey"].copy()

    merged = pd.merge(ww, ss, on="gene", how="inner")
    merged = merged.dropna(subset=["wgcna","se2","gene"])
    merged["wgcna"] = merged["wgcna"].astype(str)
    merged["se2"]   = merged["se2"].astype(str)
    merged["gene"]  = merged["gene"].astype(str)
    return merged

# ---------- Metrics ----------
def contingency_and_metrics(df: pd.DataFrame):
    M = pd.crosstab(df["wgcna"], df["se2"])  # counts
    size_w = M.sum(axis=1).to_frame("size_w")
    size_s = M.sum(axis=0).to_frame("size_s")

    J = pd.DataFrame(index=M.index, columns=M.columns, dtype=float)
    O = pd.DataFrame(index=M.index, columns=M.columns, dtype=float)
    for i in M.index:
        for j in M.columns:
            inter = M.loc[i,j]
            union = size_w.loc[i,"size_w"] + size_s.loc[j,"size_s"] - inter
            J.loc[i,j] = (inter/union) if union>0 else 0.0
            denom = min(size_w.loc[i,"size_w"], size_s.loc[j,"size_s"])
            O.loc[i,j] = (inter/denom) if denom>0 else 0.0

    # Fisher per-cell + FDR
    N = df["gene"].nunique()
    pvals = []
    for i in M.index:
        for j in M.columns:
            a = M.loc[i,j]  # intersection count
            b = size_w.loc[i,"size_w"] - a  # WGCNA module size minus intersection
            c = size_s.loc[j,"size_s"] - a  # SE2 cluster size minus intersection
            d = N - (a+b+c)  # remaining genes not in either module/cluster
            
            # Ensure all values are non-negative for Fisher's exact test
            if a < 0 or b < 0 or c < 0 or d < 0:
                print(f"Warning: Negative values in contingency table for {i}, {j}: a={a}, b={b}, c={c}, d={d}", file=sys.stderr)
                # Set p-value to 1.0 (no significance) if we can't compute it properly
                p = 1.0
            else:
                try:
                    _, p = fisher_exact([[a,b],[c,d]], alternative="greater")
                except ValueError as e:
                    print(f"Warning: Fisher's exact test failed for {i}, {j}: {e}", file=sys.stderr)
                    p = 1.0
            pvals.append(p)
    qvals = multipletests(pvals, method="fdr_bh")[1]
    Q = pd.DataFrame(qvals.reshape(len(M.index), len(M.columns)), index=M.index, columns=M.columns)
    return M, J, O, Q

def best_one_to_one(J: pd.DataFrame) -> pd.DataFrame:
    cost = 1 - J.values
    r, c = linear_sum_assignment(cost)
    return pd.DataFrame({
        "wgcna": [J.index[i] for i in r],
        "se2":   [J.columns[j] for j in c],
        "jaccard": [J.values[i,j] for i,j in zip(r,c)]
    }).sort_values("jaccard", ascending=False)

def per_module_scores(df: pd.DataFrame, pairing: pd.DataFrame) -> pd.DataFrame:
    out = []
    for _, row in pairing.iterrows():
        mw, cs = row["wgcna"], row["se2"]
        Mw = set(df.loc[df["wgcna"]==mw, "gene"])
        Cs = set(df.loc[df["se2"]==cs, "gene"])
        inter = len(Mw & Cs)
        prec = inter / len(Cs) if len(Cs)>0 else 0.0
        rec  = inter / len(Mw) if len(Mw)>0 else 0.0
        f1 = 2*prec*rec/(prec+rec) if (prec+rec)>0 else 0.0
        out.append({"wgcna":mw,"se2":cs,"precision":prec,"recall":rec,"F1":f1,"purity":rec,
                    "count_inter": inter, "size_w": len(Mw), "size_s": len(Cs)})
    return pd.DataFrame(out).sort_values("F1", ascending=False)

def pairwise_metrics(df: pd.DataFrame, max_pairs_per_module: int = 3000, seed: int = 0):
    rng = np.random.default_rng(seed)
    pos_pairs = []
    for mw, gsub in df.groupby("wgcna")["gene"]:
        gs = gsub.values
        if len(gs) >= 2:
            pairs = list(combinations(gs, 2))
            if len(pairs) > max_pairs_per_module:
                idx = rng.choice(len(pairs), size=max_pairs_per_module, replace=False)
                pairs = [pairs[i] for i in idx]
            pos_pairs.extend(pairs)
    pos = pd.DataFrame(pos_pairs, columns=["g1","g2"])
    pos["y"] = 1

    all_genes = df["gene"].unique()
    neg_pairs = set()
    target_neg = len(pos) * 2
    while len(neg_pairs) < target_neg and len(all_genes) >= 2:
        g1, g2 = rng.choice(all_genes, 2, replace=False)
        labs = df.loc[df["gene"].isin([g1,g2]), "wgcna"].unique()
        if len(labs) >= 2:
            pair = (g1,g2) if g1 < g2 else (g2,g1)
            neg_pairs.add(pair)
    neg = pd.DataFrame(list(neg_pairs), columns=["g1","g2"])
    neg["y"] = 0

    pairs = pd.concat([pos, neg], ignore_index=True)
    g2s = df.set_index("gene")["se2"].to_dict()
    pairs["yhat"] = (pairs.apply(lambda r: g2s.get(r["g1"], "_") == g2s.get(r["g2"], "_"), axis=1)).astype(int)

    mcc = matthews_corrcoef(pairs["y"], pairs["yhat"])
    score = pairs["yhat"].astype(float).values
    y = pairs["y"].values
    ap = average_precision_score(y, score)
    prec, rec, thr = precision_recall_curve(y, score)
    return {"pairs_df": pairs, "MCC": float(mcc), "AUPRC": float(ap), "prec": prec, "rec": rec}

def permutation_test(df: pd.DataFrame, n_perm: int = 50, seed: int = 0):
    rng = np.random.default_rng(seed)
    base = pairwise_metrics(df, seed=seed)
    observed_mcc = base["MCC"]
    observed_ap  = base["AUPRC"]
    null_mcc = []
    null_ap  = []
    se2_vals = df["se2"].values.copy()

    for b in range(n_perm):
        perm = se2_vals.copy()
        rng.shuffle(perm)
        df_perm = df.copy()
        df_perm["se2"] = perm
        res = pairwise_metrics(df_perm, seed=seed+b)
        null_mcc.append(res["MCC"])
        null_ap.append(res["AUPRC"])

    p_mcc = (np.sum(np.array(null_mcc) >= observed_mcc) + 1) / (n_perm + 1)
    p_ap  = (np.sum(np.array(null_ap)  >= observed_ap ) + 1) / (n_perm + 1)
    return {"observed": {"MCC": observed_mcc, "AUPRC": observed_ap},
            "null": {"MCC": null_mcc, "AUPRC": null_ap},
            "pvals": {"MCC": float(p_mcc), "AUPRC": float(p_ap)}}

# ---------- Plots ----------
def fig_heatmap_jaccard(M: pd.DataFrame, J: pd.DataFrame, Q: pd.DataFrame, out_png: str, out_pdf: str, q_thr: float = 0.05):
    fig, ax = plt.subplots(figsize=(max(6, 0.3*J.shape[1]), max(6, 0.3*J.shape[0])))
    im = ax.imshow(J.values, aspect='auto')
    ax.set_xticks(range(J.shape[1])); ax.set_xticklabels(J.columns, rotation=90)
    ax.set_yticks(range(J.shape[0])); ax.set_yticklabels(J.index)
    ax.set_xlabel("SpeakEasy2"); ax.set_ylabel("WGCNA")
    ax.set_title("Jaccard overlap (text = counts; * = q<0.05)")
    for i in range(J.shape[0]):
        for j in range(J.shape[1]):
            txt = f"{int(M.values[i,j])}"
            if Q.values[i,j] < q_thr: txt += "*"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7)
    cbar = fig.colorbar(im, ax=ax); cbar.set_label("Jaccard")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200); fig.savefig(out_pdf); plt.close(fig)

def fig_bar_per_module(scores: pd.DataFrame, out_png: str, out_pdf: str, metric: str = "purity"):
    s = scores.sort_values(metric, ascending=False)
    fig, ax = plt.subplots(figsize=(max(6, 0.25*len(s)), 6))
    ax.bar(range(len(s)), s[metric].values)
    ax.set_xticks(range(len(s))); ax.set_xticklabels(s["wgcna"].values, rotation=90)
    ax.set_ylabel(metric); ax.set_title(f"{metric} per WGCNA module (matched by Hungarian)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200); fig.savefig(out_pdf); plt.close(fig)

def fig_pr_curve(res_pairs: Dict, out_png: str, out_pdf: str):
    prec = res_pairs["prec"]; rec  = res_pairs["rec"]; ap   = res_pairs["AUPRC"]
    fig, ax = plt.subplots(figsize=(6,5))
    ax.plot(rec, prec)
    ax.set_xlabel("Recall"); ax.set_ylabel("Precision"); ax.set_title(f"Precision-Recall (AUPRC={ap:.3f})")
    fig.tight_layout(); fig.savefig(out_png, dpi=200); fig.savefig(out_pdf); plt.close(fig)

def fig_null_hist(observed: float, null_vals, label: str, out_png: str, out_pdf: str):
    fig, ax = plt.subplots(figsize=(6,4))
    ax.hist(null_vals, bins=30)
    ax.axvline(observed, linestyle="--")
    ax.set_title(f"Null distribution for {label} (permute SE2 labels)")
    ax.set_xlabel(label); ax.set_ylabel("Count")
    fig.tight_layout(); fig.savefig(out_png, dpi=200); fig.savefig(out_pdf); plt.close(fig)

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wgcna", required=True)
    ap.add_argument("--se2", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--gene-col-wgcna", default=None)
    ap.add_argument("--module-col", default=None)
    ap.add_argument("--gene-col-se2", default=None)
    ap.add_argument("--se2-col", default=None)
    ap.add_argument("--drop-grey", action="store_true")
    ap.add_argument("--max-pairs-per-module", type=int, default=3000)
    ap.add_argument("--n-perm", type=int, default=50)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    w = smart_read(args.wgcna)
    s = smart_read(args.se2)

    try:
        df = standardize(w, s,
                         gene_col_w=args.gene_col_wgcna,
                         wgcna_col=args.module_col,
                         gene_col_s=args.gene_col_se2,
                         se2_col=args.se2_col,
                         drop_grey=args.drop_grey)
    except Exception as e:
        with open(os.path.join(args.outdir, "ERROR.txt"), "w", encoding="utf-8") as f:
            f.write(str(e))
        print("ERROR:", e, file=sys.stderr)
        sys.exit(1)

    # Save clean labels
    df.to_csv(os.path.join(args.outdir, "clean_labels.csv"), index=False)

    # Core matrices
    M, J, O, Q = contingency_and_metrics(df)
    M.to_csv(os.path.join(args.outdir, "contingency_counts.csv"))
    J.to_csv(os.path.join(args.outdir, "jaccard.csv"))
    O.to_csv(os.path.join(args.outdir, "overlap.csv"))
    Q.to_csv(os.path.join(args.outdir, "qvalues_fdr.csv"))

    # Hungarian and per-module scores
    pairs = best_one_to_one(J)
    pairs.to_csv(os.path.join(args.outdir, "best_pairs_hungarian.csv"), index=False)
    scores = per_module_scores(df, pairs)
    scores.to_csv(os.path.join(args.outdir, "per_module_scores.csv"), index=False)

    # Pairwise and permutations
    res_pairs = pairwise_metrics(df, max_pairs_per_module=args.max_pairs_per_module, seed=0)
    perm = permutation_test(df, n_perm=args.n_perm, seed=0)
    pd.DataFrame({"MCC": perm["null"]["MCC"]}).to_csv(os.path.join(args.outdir,"perm_null_MCC.csv"), index=False)
    pd.DataFrame({"AUPRC": perm["null"]["AUPRC"]}).to_csv(os.path.join(args.outdir,"perm_null_AUPRC.csv"), index=False)

    summary = {
        "n_genes": int(df["gene"].nunique()),
        "n_wgcna_modules": int(df["wgcna"].nunique()),
        "n_se2_clusters": int(df["se2"].nunique()),
        "hungarian_mean_jaccard": float(pairs["jaccard"].mean()) if len(pairs)>0 else None,
        "mean_module_purity": float(scores["purity"].mean()) if len(scores)>0 else None,
        "pairwise_MCC": float(res_pairs["MCC"]),
        "pairwise_AUPRC": float(res_pairs["AUPRC"]),
        "perm_pvals": perm["pvals"]
    }
    with open(os.path.join(args.outdir, "summary_metrics.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # Plots
    fig_heatmap_jaccard(M, J, Q,
        out_png=os.path.join(args.outdir, "heatmap_jaccard.png"),
        out_pdf=os.path.join(args.outdir, "heatmap_jaccard.pdf"))
    fig_bar_per_module(scores,
        out_png=os.path.join(args.outdir, "per_module_purity.png"),
        out_pdf=os.path.join(args.outdir, "per_module_purity.pdf"), metric="purity")
    fig_bar_per_module(scores,
        out_png=os.path.join(args.outdir, "per_module_F1.png"),
        out_pdf=os.path.join(args.outdir, "per_module_F1.pdf"), metric="F1")
    fig_pr_curve(res_pairs,
        out_png=os.path.join(args.outdir, "precision_recall_curve.png"),
        out_pdf=os.path.join(args.outdir, "precision_recall_curve.pdf"))
    fig_null_hist(perm["observed"]["MCC"], perm["null"]["MCC"], "MCC",
        out_png=os.path.join(args.outdir, "null_hist_MCC.png"),
        out_pdf=os.path.join(args.outdir, "null_hist_MCC.pdf"))
    fig_null_hist(perm["observed"]["AUPRC"], perm["null"]["AUPRC"], "AUPRC",
        out_png=os.path.join(args.outdir, "null_hist_AUPRC.png"),
        out_pdf=os.path.join(args.outdir, "null_hist_AUPRC.pdf"))

    # README
    readme = """# Outputs: WGCNA vs SpeakEasy2 overlap

**Clean labels**
- `clean_labels.csv` — טבלת בסיס עם `gene,wgcna,se2` אחרי איחוד וניקוי.

**חפיפות פר-מודול**
- `contingency_counts.csv` — מטריצת ספירות |Mi ∩ Sj|.
- `jaccard.csv` — מטריצת Jaccard בין כל מודול WGCNA לכל קהילת SE2.
- `overlap.csv` — מטריצת Overlap (Szymkiewicz–Simpson).
- `qvalues_fdr.csv` — ערכי q (FDR) למבחני Fisher per-cell.

**התאמה אחד-לאחד**
- `best_pairs_hungarian.csv` — צמדי מודול↔קהילה שנבחרו ע"י Hungarian עם ערך Jaccard.

**ביצועי מודולים מול הקהילה התואמת**
- `per_module_scores.csv` — precision/recall/F1/purity + גדלים וחיתוך.

**מדדי זוג-גנים (הוכחה כללית)**
- `summary_metrics.json` — סיכום כולל MCC/AUPRC ו-p-values מהפרמוטציות.
- `perm_null_MCC.csv`, `perm_null_AUPRC.csv` — התפלגויות null.

**ויזואליזציות**
- `heatmap_jaccard.(png|pdf)` — Heatmap Jaccard (טקסט: ספירות; *: q<0.05).
- `per_module_purity.(png|pdf)` — Purity לפר מודול.
- `per_module_F1.(png|pdf)` — F1 לפר מודול.
- `precision_recall_curve.(png|pdf)` — עקומת Precision–Recall (AUPRC).
- `null_hist_MCC.(png|pdf)`, `null_hist_AUPRC.(png|pdf)` — היסטוגרמות null מול הערך הנצפה.

## פרשנות מהירה
- אלכסון "חם" ב-Heatmap + q<0.05 רבים ⇒ חפיפה עקבית.
- Purity/Recall/F1 גבוהים ⇒ מודול WGCNA מכוסה היטב בקהילת SE2 מתאימה.
- MCC/AUPRC גבוהים משמעותית מ-null ⇒ שני הגדרות הקלאסטרים נותנות תמונה קוהרנטית ברמת זוג-גנים.
"""
    with open(os.path.join(args.outdir, "README.md"), "w", encoding="utf-8") as f:
        f.write(readme)

    print("OK")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
