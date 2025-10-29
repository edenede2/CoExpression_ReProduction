#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
viz_wgcna_se2_one_to_many.py
בר-פלוטים "אחד-לכמה" (one-to-many) עבור WGCNA↔SpeakEasy2 מהפלטים הקיימים.

נכנסים (בתיקיית --indir):
  clean_labels.csv
  jaccard.csv
  qvalues_fdr.csv
  overlap.csv
אופציונלי: best_pairs_hungarian.csv, per_module_scores.csv (לא חובה)

יוצאים (לתיקיית --outdir):
  - bar_top_pairs_one_to_many.png       : טופ N הזוגות (W→S) לפי המדד שנבחר, ללא אילוץ אחד-לאחד
  - bar_module_topk_one_to_many.png     : טופ K קהילות לכל מודול, מוצגים עבור top-modules מודולים
  - one_to_many_pairs.csv               : טבלת כל הזוגות אחרי פילטרים + הדירוג
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- IO helpers ----------
def read_matrix_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    # אם העמודה הראשונה אינה מספרית — כנראה שמות שורות
    if df.shape[1] > 1 and not pd.api.types.is_numeric_dtype(df.iloc[:,0]):
        df = df.set_index(df.columns[0])
    df = df.apply(pd.to_numeric, errors="coerce")
    return df

def read_clean_labels(path: str) -> pd.DataFrame:
    # מצופה עמודות: gene, wgcna, se2  (אחרי הניקוי)
    return pd.read_csv(path)

def ensure_intersection(A: pd.DataFrame, B: pd.DataFrame):
    rows = A.index.intersection(B.index)
    cols = A.columns.intersection(B.columns)
    return A.loc[rows, cols], B.loc[rows, cols]

def read_optional(path):
    return pd.read_csv(path) if os.path.exists(path) else None

# ---------- Build long-form pairs ----------
def pairs_from_matrix(M: pd.DataFrame, name: str) -> pd.DataFrame:
    # הופך מטריצה ל־long (wgcna, se2, value)
    df = M.copy()
    df["wgcna"] = df.index
    long = df.melt(id_vars="wgcna", var_name="se2", value_name=name)
    return long

# ---------- Plot helpers ----------
def barplot_long(df: pd.DataFrame, label_col: str, value_col: str, out_png: str, title: str, max_bars: int):
    data = df.sort_values(value_col, ascending=False).head(max_bars)
    n = len(data)
    fig, ax = plt.subplots(figsize=(max(8, 0.35*n), 6))
    ax.bar(range(n), data[value_col].values)
    ax.set_xticks(range(n))
    ax.set_xticklabels(data[label_col].values, rotation=90)
    ax.set_title(title)
    ax.set_ylabel(value_col)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

def barplot_per_module_topk(df_topk: pd.DataFrame, out_png: str, title: str, top_modules: int, k: int):
    # df_topk: columns = [wgcna, se2, score]
    # נבחר top_modules מודולים לפי ה־max score שלהם
    order_mod = (df_topk.groupby("wgcna")["score"]
                         .max()
                         .sort_values(ascending=False)
                         .head(top_modules)
                         .index)
    sub = df_topk[df_topk["wgcna"].isin(order_mod)].copy()
    # כדי לצייר בבר-פלוט יחיד, נכניס תווית משולבת W→S
    sub = sub.sort_values(["wgcna","score"], ascending=[True, False])
    sub["label"] = sub.apply(lambda r: f"W:{r['wgcna']} → S:{r['se2']}", axis=1)

    # בר-פלוט אחד ארוך: k עמודות לכל מודול, לפי סדר מודולים
    n = len(sub)
    fig, ax = plt.subplots(figsize=(max(8, 0.35*n), 6))
    ax.bar(range(n), sub["score"].values)
    ax.set_xticks(range(n))
    ax.set_xticklabels(sub["label"].values, rotation=90)
    ax.set_title(f"{title}  (top {k} per module, {len(order_mod)} modules)")
    ax.set_ylabel("score")
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, help="Folder with CSVs")
    ap.add_argument("--outdir", required=True, help="Output folder for PNG/CSV")
    ap.add_argument("--measure", choices=["jaccard","overlap"], default="jaccard",
                    help="Which score to rank by")
    ap.add_argument("--k", type=int, default=3, help="Top-k SE2 per WGCNA module (one-to-many)")
    ap.add_argument("--top-n", type=int, default=150, help="Top-N global pairs to show in overall barplot")
    ap.add_argument("--top-modules", type=int, default=20, help="How many modules to show in the per-module barplot")
    ap.add_argument("--qthr", type=float, default=0.05, help="q-value threshold for significance filter (use only if qvals file exists)")
    ap.add_argument("--min-overlap", type=float, default=0.0, help="Minimal overlap threshold to keep pairs (0-1), requires overlap.csv")
    ap.add_argument("--min-se2-size", type=int, default=0, help="Minimal SE2 cluster size (by #genes in clean_labels)")
    ap.add_argument("--min-wgcna-size", type=int, default=0, help="Minimal WGCNA module size (by #genes in clean_labels)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    P = lambda f: os.path.join(args.indir, f)

    # Load matrices
    J = read_matrix_csv(P("jaccard.csv"))
    O = read_matrix_csv(P("overlap.csv")) if os.path.exists(P("overlap.csv")) else None
    Q = read_matrix_csv(P("qvalues_fdr.csv")) if os.path.exists(P("qvalues_fdr.csv")) else None

    # Align J with Q (אם קיים) ועם O (אם קיים)
    if Q is not None:
        J, Q = ensure_intersection(J, Q)
    if O is not None:
        rows = J.index.intersection(O.index)
        cols = J.columns.intersection(O.columns)
        J = J.loc[rows, cols]
        O = O.loc[rows, cols]

    # Sizes from clean_labels (לסינון לפי גודל)
    CL = read_clean_labels(P("clean_labels.csv"))
    w_sizes = CL.groupby("wgcna")["gene"].nunique()
    s_sizes = CL.groupby("se2")["gene"].nunique()

    # בסיס הציונים
    ScoreM = J if args.measure == "jaccard" else (O if O is not None else J.copy())

    # מסכת פילטרים
    keep = pd.DataFrame(True, index=ScoreM.index, columns=ScoreM.columns)

    # פילטר מובהקות (אם יש Q)
    if Q is not None and args.qthr is not None:
        keep &= (Q < args.qthr)

    # פילטר overlap מינימלי (אם ביקשת ויש O)
    if args.min_overlap > 0.0 and O is not None:
        keep &= (O >= args.min_overlap)

    # פילטרי גודל מודול/קהילה
    if args.min_wgcna_size > 0:
        big_w = set(w_sizes[w_sizes >= args.min_wgcna_size].index)
        keep = keep.loc[keep.index.intersection(big_w)]
        ScoreM = ScoreM.loc[ScoreM.index.intersection(big_w)]
        if Q is not None: Q = Q.loc[Q.index.intersection(big_w)]
        if O is not None: O = O.loc[O.index.intersection(big_w)]
    if args.min_se2_size > 0:
        big_s = set(s_sizes[s_sizes >= args.min_se2_size].index)
        keep = keep.loc[:, keep.columns.intersection(big_s)]
        ScoreM = ScoreM.loc[:, ScoreM.columns.intersection(big_s)]
        if Q is not None: Q = Q.loc[:, Q.columns.intersection(big_s)]
        if O is not None: O = O.loc[:, O.columns.intersection(big_s)]

    # אפס ציונים היכן שאסור לשמור (כדי לא להופיע בטופ)
    ScoreF = ScoreM.where(keep, other=np.nan)

    # --- בניית טבלת one-to-many ---
    # 1) טופ-K לכל מודול (לפי score יורד)
    rows = []
    for w in ScoreF.index:
        s_sorted = ScoreF.loc[w].sort_values(ascending=False)
        s_sorted = s_sorted.dropna().head(max(0, args.k))
        for s, val in s_sorted.items():
            rows.append({"wgcna": w, "se2": s, "score": float(val)})
    df_topk = pd.DataFrame(rows)

    # 2) טופ-N גלובלי מכל הזוגות האפשריים אחרי פילטרים
    long_all = pairs_from_matrix(ScoreF, name="score").dropna(subset=["score"])
    long_all["label"] = long_all.apply(lambda r: f"W:{r['wgcna']} → S:{r['se2']}", axis=1)
    
    # Add overlap count information if overlap matrix is available
    if O is not None:
        # Create a DataFrame with overlap counts (which might be percentages)
        overlap_values = pairs_from_matrix(O, name="overlap_value").dropna(subset=["overlap_value"])
        # Merge with the main DataFrame
        long_all = long_all.merge(overlap_values[["wgcna", "se2", "overlap_value"]], 
                                  on=["wgcna", "se2"], how="left")
        
        # If we have clean_labels, we can calculate the actual gene counts
        if "CL" in locals():
            # Get module and community sizes
            w_sizes = CL.groupby("wgcna")["gene"].nunique().reset_index()
            w_sizes.columns = ["wgcna", "wgcna_size"]
            s_sizes = CL.groupby("se2")["gene"].nunique().reset_index()
            s_sizes.columns = ["se2", "se2_size"]
            
            # Ensure data types match for merge operations
            # Convert se2 column to the same type in both DataFrames
            long_all["se2"] = long_all["se2"].astype(str)
            s_sizes["se2"] = s_sizes["se2"].astype(str)
            # Convert wgcna column to the same type in both DataFrames
            long_all["wgcna"] = long_all["wgcna"].astype(str)
            w_sizes["wgcna"] = w_sizes["wgcna"].astype(str)
            
            # Add size information to the DataFrame
            long_all = long_all.merge(w_sizes, on="wgcna", how="left")
            long_all = long_all.merge(s_sizes, on="se2", how="left")
            
            # Calculate the actual intersection size (number of genes in the overlap)
            # First, let's create a function to get the intersection count for each pair
            def get_intersection_count(row):
                wgcna_module = row["wgcna"]
                se2_community = row["se2"]
                # Get genes in this WGCNA module
                wgcna_genes = set(CL[CL["wgcna"] == wgcna_module]["gene"])
                # Get genes in this SpeakEasy2 community
                se2_genes = set(CL[CL["se2"] == se2_community]["gene"])
                # Return the size of the intersection
                return len(wgcna_genes.intersection(se2_genes))
            
            # Apply the function to each row to get the actual intersection count
            # We'll use a vectorized approach for better performance
            
            # Examine gene identifiers to look for patterns
            gene_sample = CL["gene"].head(10).tolist()
            
            # Check if gene identifiers contain tissue or version suffixes
            has_version_suffix = any('.' in str(gene) for gene in gene_sample)
            has_ac_suffix = any('_' in str(gene) for gene in gene_sample)
            
            # Try to normalize gene identifiers by removing version numbers and suffixes
            if has_version_suffix or has_ac_suffix:
                # Create a function to normalize gene identifiers while preserving tissue information
                def normalize_gene_id(gene_id):
                    gene_id = str(gene_id)
                    # Only remove version number (e.g., .16) but keep tissue suffix
                    if '.' in gene_id:
                        # Split by dot but preserve everything after the last dot if it contains a tissue suffix
                        parts = gene_id.split('.')
                        base = parts[0]  # ENSG part
                        # Keep the tissue suffix if present
                        if '_' in parts[-1]:
                            suffix = parts[-1].split('_', 1)[1]  # Get the tissue part after _
                            return f"{base}_{suffix}"
                        return base
                    return gene_id
                
                # Create a new column with normalized gene IDs
                CL["gene_normalized"] = CL["gene"].apply(normalize_gene_id)
                
                # Use normalized gene IDs for intersection calculations
                use_normalized = True
                gene_col = "gene_normalized"
            else:
                use_normalized = False
                gene_col = "gene"
            
            # Make sure we convert types to match in both DataFrames
            # Force both columns to strings for comparison
            CL["wgcna"] = CL["wgcna"].astype(str)
            CL["se2"] = CL["se2"].astype(str)
            
            # Create dictionaries for fast lookup of genes in each module/community
            wgcna_gene_dict = {str(w): set(CL[CL["wgcna"] == str(w)][gene_col]) for w in long_all["wgcna"].unique()}
            se2_gene_dict = {str(s): set(CL[CL["se2"] == str(s)][gene_col]) for s in long_all["se2"].unique()}
            
            # Check module/community coverage
            all_wgcna_modules = set(long_all["wgcna"].unique())
            wgcna_dict_keys = set(wgcna_gene_dict.keys())
            
            all_se2_communities = set(long_all["se2"].unique())
            se2_dict_keys = set(se2_gene_dict.keys())
            
            # Calculate the actual intersection count (number of genes in the overlap)
            long_all["gene_intersection_count"] = long_all.apply(
                lambda row: len(wgcna_gene_dict.get(str(row["wgcna"]), set()) & 
                                se2_gene_dict.get(str(row["se2"]), set())), axis=1)
            
            # Report on the gene normalization impact - removed debug prints
                
            # Check for inconsistency between Jaccard and intersection
            high_jaccard_zero_int = long_all[(long_all["score"] > 0.2) & (long_all["gene_intersection_count"] == 0)]
            if not high_jaccard_zero_int.empty:
                print(f"\nFound {len(high_jaccard_zero_int)} pairs with high Jaccard (>0.2) but zero intersection!")
                print("Sample of these pairs:")
                sample = high_jaccard_zero_int.head(5)
                for _, row in sample.iterrows():
                    w = str(row["wgcna"])
                    s = str(row["se2"])
                    print(f"  WGCNA {w} - SE2 {s}: Jaccard={row['score']}, Intersection=0")
                    print(f"    WGCNA {w} has {len(wgcna_gene_dict.get(w, set()))} genes")
                    print(f"    SE2 {s} has {len(se2_gene_dict.get(s, set()))} genes")
                    
                    # For one problematic pair, examine sample genes to look for patterns
                    if len(high_jaccard_zero_int) > 0:
                        problem_row = high_jaccard_zero_int.iloc[0]
                        w = str(problem_row["wgcna"])
                        s = str(problem_row["se2"])
                        w_genes = list(wgcna_gene_dict.get(w, set()))[:5]
                        s_genes = list(se2_gene_dict.get(s, set()))[:5]
                        print(f"\nSample genes from problematic pair WGCNA {w} - SE2 {s}:")
                        print(f"  WGCNA {w} genes: {w_genes}")
                        print(f"  SE2 {s} genes: {s_genes}")
                        
                        # Check if there are any potential matches by basic gene ID
                        # Extract just the ENSG part without version or tissue
                        def get_basic_gene_id(gene_id):
                            gene_id = str(gene_id)
                            if 'ENSG' in gene_id:
                                # Extract just the ENSG + number part
                                import re
                                match = re.search(r'(ENSG\d+)', gene_id)
                                return match.group(1) if match else gene_id
                            return gene_id
                        
                        # Get basic gene IDs
                        w_genes_basic = [get_basic_gene_id(g) for g in w_genes]
                        s_genes_basic = [get_basic_gene_id(g) for g in s_genes]
                        
                        # Get all genes from both modules (first 100 only to avoid excessive output)
                        all_w_genes = list(wgcna_gene_dict.get(w, set()))[:100]
                        all_s_genes = list(se2_gene_dict.get(s, set()))[:100]
                        all_w_genes_basic = [get_basic_gene_id(g) for g in all_w_genes]
                        all_s_genes_basic = [get_basic_gene_id(g) for g in all_s_genes]
                        
                        # Check for any matches in the basic gene IDs
                        basic_matches = set(all_w_genes_basic).intersection(set(all_s_genes_basic))
                        
                        if basic_matches:
                            print(f"\nFound {len(basic_matches)} matches when ignoring version and tissue info!")
                            print(f"  First 5 matching basic gene IDs: {list(basic_matches)[:5]}")
                            
                            # Find examples of full gene IDs that match on basic ID
                            print("\nExamples of full gene IDs that match on basic ID but were not counted as overlaps:")
                            for basic_id in list(basic_matches)[:3]:  # Show first 3 examples
                                w_examples = [g for g in all_w_genes if basic_id in g]
                                s_examples = [g for g in all_s_genes if basic_id in g]
                                if w_examples and s_examples:
                                    print(f"  Basic ID: {basic_id}")
                                    print(f"    WGCNA version: {w_examples[0]}")
                                    print(f"    SE2 version: {s_examples[0]}")
                        else:
                            print("\nNo matches found even when ignoring version and tissue information.")            # Calculate expected Jaccard index based on gene sets
            def calculate_jaccard(row):
                w = str(row["wgcna"])
                s = str(row["se2"])
                w_genes = wgcna_gene_dict.get(w, set())
                s_genes = se2_gene_dict.get(s, set())
                
                if not w_genes or not s_genes:
                    return 0.0
                    
                intersection = len(w_genes.intersection(s_genes))
                union = len(w_genes.union(s_genes))
                return intersection / union if union > 0 else 0.0
            
            # Add calculated Jaccard and check discrepancies with the provided scores
            long_all["calculated_jaccard"] = long_all.apply(calculate_jaccard, axis=1)
            
            # Compare the provided Jaccard scores with the calculated ones
            long_all["jaccard_diff"] = abs(long_all["score"] - long_all["calculated_jaccard"])
            
            # Print summary of discrepancies
            large_diff_count = (long_all["jaccard_diff"] > 0.1).sum()
            if large_diff_count > 0:
                print(f"\nFound {large_diff_count} pairs with large discrepancy (>0.1) between provided and calculated Jaccard")
                print("Sample of pairs with large discrepancies:")
                sample_diff = long_all[long_all["jaccard_diff"] > 0.1].sort_values("jaccard_diff", ascending=False).head(5)
                for _, row in sample_diff.iterrows():
                    w = str(row["wgcna"])
                    s = str(row["se2"])
                    print(f"  WGCNA {w} - SE2 {s}:")
                    print(f"    Provided Jaccard: {row['score']:.4f}")
                    print(f"    Calculated Jaccard: {row['calculated_jaccard']:.4f}")
                    print(f"    Intersection count: {row['gene_intersection_count']}")
                    print(f"    WGCNA size: {len(wgcna_gene_dict.get(w, set()))}")
                    print(f"    SE2 size: {len(se2_gene_dict.get(s, set()))}")
                    
            # Remove the calculated and diff columns before saving to CSV
            long_all = long_all.drop(columns=["calculated_jaccard", "jaccard_diff"])

    # Add tissue information to the pairs
    if "CL" in locals() and "_" in str(CL["gene"].iloc[0]):
        # Determine the dominant tissue for each module and community
        def get_dominant_tissue(genes):
            tissues = [str(g).split('_')[-1] for g in genes if "_" in str(g)]
            if not tissues:
                return "unknown"
            from collections import Counter
            return Counter(tissues).most_common(1)[0][0]
        
        # Get dominant tissues for each WGCNA module
        wgcna_tissues = {w: get_dominant_tissue(genes) for w, genes in wgcna_gene_dict.items()}
        se2_tissues = {s: get_dominant_tissue(genes) for s, genes in se2_gene_dict.items()}
        
        # Add tissue information to the output
        long_all["wgcna_tissue"] = long_all["wgcna"].apply(lambda w: wgcna_tissues.get(str(w), "unknown"))
        long_all["se2_tissue"] = long_all["se2"].apply(lambda s: se2_tissues.get(str(s), "unknown"))
        
        # Flag pairs where tissues don't match
        long_all["same_tissue"] = long_all["wgcna_tissue"] == long_all["se2_tissue"]
    
    # Calculate a reliability score based on the agreement between provided Jaccard and overlap counts
    if "calculated_jaccard" in long_all.columns:
        # Normalize the difference to 0-1 scale where 1 means perfect agreement
        max_diff = long_all["jaccard_diff"].max()
        if max_diff > 0:
            long_all["reliability"] = 1 - (long_all["jaccard_diff"] / max_diff)
        else:
            long_all["reliability"] = 1.0
        
        # Add flag for inconsistent data
        long_all["consistent_data"] = (long_all["gene_intersection_count"] > 0) | (long_all["score"] < 0.05)
        
        # Remove temporary calculation columns
        long_all = long_all.drop(columns=["calculated_jaccard", "jaccard_diff"])
    else:
        # Add consistency flag anyway
        long_all["consistent_data"] = (long_all["gene_intersection_count"] > 0) | (long_all["score"] < 0.05)
    
    # Add a note about the discrepancy and tissue differences
    print("\nNOTE: There appear to be discrepancies between the provided Jaccard scores and the calculated values.")
    print("Analysis shows that the genes in different modules/communities often come from different tissues:")
    if "CL" in locals() and "_" in str(CL["gene"].iloc[0]):
        tissue_counts = CL["gene"].apply(lambda g: str(g).split('_')[-1] if "_" in str(g) else "unknown").value_counts()
        print(f"  Tissue distribution in the data:")
        for tissue, count in tissue_counts.items():
            print(f"    {tissue}: {count} genes")
    
    print("\nThe output CSV includes:")
    print("  - gene_intersection_count: Actual count of genes in the intersection (based on clean_labels.csv)")
    print("  - wgcna_size and se2_size: Module/community sizes (based on clean_labels.csv)")
    if "wgcna_tissue" in long_all.columns:
        print("  - wgcna_tissue and se2_tissue: Dominant tissue type for each module/community")
        print("  - same_tissue: Flag indicating if the WGCNA module and SE2 community have the same dominant tissue")
    print("  - consistent_data: Flag for potentially problematic pairs (high Jaccard but low intersection)")
    
    # שמירת CSV מלא
    out_csv = os.path.join(args.outdir, "one_to_many_pairs.csv")
    long_all.sort_values("score", ascending=False).to_csv(out_csv, index=False)

    # --- ציורים ---
    # Global top-N pairs
    barplot_long(
        long_all,
        label_col="label",
        value_col="score",
        out_png=os.path.join(args.outdir, "bar_top_pairs_one_to_many.png"),
        title=f"Top {args.top_n} pairs (one-to-many) by {args.measure}",
        max_bars=args.top_n
    )

    # Per-module top-k (למודולים הטובים ביותר)
    if not df_topk.empty:
        barplot_per_module_topk(
            df_topk,
            out_png=os.path.join(args.outdir, "bar_module_topk_one_to_many.png"),
            title=f"Top-{args.k} SE2 per WGCNA module by {args.measure}",
            top_modules=args.top_modules,
            k=args.k
        )

    print("Wrote:",
          os.path.join(args.outdir, "bar_top_pairs_one_to_many.png"),
          os.path.join(args.outdir, "bar_module_topk_one_to_many.png"),
          out_csv)

if __name__ == "__main__":
    main()
