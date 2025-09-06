import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
import seaborn as sns
import matplotlib.patheffects as pe
from mpl_toolkits.axes_grid1.inset_locator import inset_axes



def add_bottom_module_strip(ax_hm, module_ids, fontsize=6, rotation=90, frac=0.16):

    strip = inset_axes(ax_hm, width="100%", height=f"{frac*100}%",
                       loc="lower center",
                       bbox_to_anchor=(0, -frac, 1, frac),
                       bbox_transform=ax_hm.transAxes,
                       borderpad=0)
    n = len(module_ids)
    strip.set_xlim(-0.5, n - 0.5)
    strip.set_ylim(0, 1)
    strip.axis("off")

    x = np.arange(n)
    for j, mid in enumerate(module_ids):
        strip.text(j, 0.0, str(int(mid)), ha="center", va="bottom",
                   rotation=rotation, fontsize=fontsize, clip_on=False)
    return strip



def _std_cols(df):
    low = {c.lower(): c for c in df.columns}
    def pick(*keys):
        for k in keys:
            if k in low: return low[k]
        return None
    return {
        "Cluster":    pick("cluster"),
        "Description":pick("description","pathway","term","name"),
        "ID":         pick("id"),
        "category":   pick("category","catagory"),
        "subcategory":pick("subcategory","subcatagory","sub_category","sub-cat")
    }

def _extract_row_meta_from_kegg(kegg_long_csv, row_labels):
    df = pd.read_csv(kegg_long_csv, engine="python", on_bad_lines="skip")
    cols = _std_cols(df)
    if cols["Description"] is None:
        return pd.DataFrame(index=row_labels, columns=["category","subcategory"])

    desc_col = cols["Description"]; cat_col = cols["category"]; sub_col = cols["subcategory"]

    meta = (df.dropna(subset=[desc_col])
              .drop_duplicates([desc_col])
              .rename(columns={
                  desc_col:"Description",
                  (cat_col or "category"):"category",
                  (sub_col or "subcategory"):"subcategory"
              })[["Description","category","subcategory"]])

    def parse_desc(lbl):
        if "—" in lbl: return lbl.split("—", 1)[-1].strip()
        if " - " in lbl: return lbl.split(" - ", 1)[-1].strip()
        return lbl.strip()

    desc_series = pd.Series([parse_desc(x) for x in row_labels], index=row_labels, name="Description")
    row_meta = desc_series.to_frame().merge(meta, on="Description", how="left").set_index(desc_series.index)
    return row_meta[["category","subcategory"]]

def build_annotation_palettes_from_kegg(kegg_long_csv, base_pal=None):
    df = pd.read_csv(kegg_long_csv, engine="python", on_bad_lines="skip")
    cols = _std_cols(df)
    cat = df[cols["category"]] if cols["category"] in df.columns else pd.Series(dtype=object)
    sub = df[cols["subcategory"]] if cols["subcategory"] in df.columns else pd.Series(dtype=object)

    def _uniq(vals):
        return [x for x in pd.unique(vals.astype(str).str.strip()) if (x!='nan') and (x!='None')]

    uniq_cat = _uniq(cat) if not cat.empty else []
    uniq_sub = _uniq(sub) if not sub.empty else []

    if base_pal is None:
        pool = (list(plt.cm.tab20.colors)+list(plt.cm.tab20b.colors)+
                list(plt.cm.tab20c.colors)+list(plt.cm.Set3.colors))
    else:
        pool = list(base_pal)

    palette_cat = {lab: pool[i % len(pool)] for i, lab in enumerate(uniq_cat)}
    palette_sub = {lab: pool[i % len(pool)] for i, lab in enumerate(uniq_sub)}
    return palette_cat, palette_sub

def _cluster_orders(logp_df, row_meta,
                    cluster_cols=True, col_metric="euclidean", col_method="average",
                    cluster_rows=True, row_metric="euclidean", row_method="average",
                    row_grouping=("category","subcategory")):
    if cluster_cols:
        gcol = sns.clustermap(
            np.nan_to_num(logp_df.values, nan=0.0),
            row_cluster=False, col_cluster=True,
            metric=col_metric, method=col_method,
            cmap="Blues", cbar=False, xticklabels=False, yticklabels=False
        )
        col_order = [logp_df.columns[i] for i in gcol.dendrogram_col.reordered_ind]
        plt.close(gcol.fig)
    else:
        col_order = list(logp_df.columns)

    rm = row_meta.copy()
    def _norm(x): return np.nan if pd.isna(x) else str(x).strip()
    rm["category"]    = rm["category"].map(_norm)
    rm["subcategory"] = rm["subcategory"].map(_norm)

    sep_positions = [] 

    if not cluster_rows:
        row_order = list(logp_df.index)
    else:
        if row_grouping is None:
            grow = sns.clustermap(
                np.nan_to_num(logp_df.values, nan=0.0),
                row_cluster=True, col_cluster=False,
                metric=row_metric, method=row_method,
                cmap="Blues", cbar=False, xticklabels=False, yticklabels=False
            )
            row_order = [logp_df.index[i] for i in grow.dendrogram_row.reordered_ind]
            plt.close(grow.fig)
        else:
            keys = list(row_grouping)
            for k in keys:
                cats = pd.unique(rm[k].dropna())
                rm[f"__ord_{k}__"] = pd.Categorical(rm[k], categories=cats, ordered=True)

            sorted_rm = rm.sort_values(by=[f"__ord_{k}__" for k in keys],
                                       na_position="last", kind="mergesort")

            row_order = []
            last_cat = None; cum_len = 0
            for grp_keys, sub_rm in sorted_rm.groupby(keys, sort=False):
                idxs = sub_rm.index
                sub  = logp_df.loc[idxs]
                if len(sub) > 1:
                    g = sns.clustermap(
                        np.nan_to_num(sub.values, nan=0.0),
                        row_cluster=True, col_cluster=False,
                        metric=row_metric, method=row_method,
                        cmap="Blues", cbar=False, xticklabels=False, yticklabels=False
                    )
                    reord = [sub.index[i] for i in g.dendrogram_row.reordered_ind]
                    plt.close(g.fig)
                else:
                    reord = list(idxs)
                curr_cat = sub_rm["category"].iloc[0]
                if last_cat is None:
                    last_cat = curr_cat
                elif curr_cat != last_cat:
                    sep_positions.append(cum_len - 0.5) 
                    last_cat = curr_cat
                row_order.extend(reord)
                cum_len += len(reord)

    return row_order, col_order, sep_positions


def _read_module_sizes(modules_details_tsv):

    md = pd.read_csv(modules_details_tsv, sep="\t")
    id_col = "Cluster ID"
    size_col = None
    for cand in ["Cluster Size","Size","Genes","Gene Count","N_genes","total"]:
        if cand in md.columns:
            size_col = cand; break

    if size_col is None and all(c in md.columns for c in ["AC","MF","PCG"]):
        sizes = (pd.to_numeric(md["AC"], errors="coerce").fillna(0) +
                 pd.to_numeric(md["MF"], errors="coerce").fillna(0) +
                 pd.to_numeric(md["PCG"], errors="coerce").fillna(0))
    elif size_col is not None:
        sizes = pd.to_numeric(md[size_col], errors="coerce").fillna(0)
    else:
        sizes = pd.Series([0]*len(md))

    out = pd.Series(sizes.values, index=pd.to_numeric(md[id_col], errors="coerce").astype("Int64")).dropna()
    out.index = out.index.astype(int)
    return out

def export_kegg_ts_ct_enrichment_pdf(
    kegg_long_csv,
    modules_details_tsv,
    pdf_path="kegg_TS_CT_enrichment.pdf",
    selection="global_top", K_GLOBAL=60, N_PER_MODULE=4, P_THRESH=0.05,
    cap=6.0, star_from="p",
    cluster_cols=True,  col_metric="euclidean", col_method="average",
    cluster_rows=True,  row_metric="euclidean", row_method="average",
    row_grouping=("category","subcategory"),
    show_row_labels=False, legend_fontsize=8,
    star_color="crimson", star_fontsize=10, star_outline=True,
    strip_width=0.018,
    add_hist=True, bar_color="#4477AA", bar_width=0.85, hist_height_ratio=0.18,
    hist_stacked_regions=False,
    show_region_legend=None,        
    hist_regions_cols=("AC","MF","PCG"),
    hist_region_colors=None,   
    palette_cat=None, palette_sub=None,
    build_matrix_fn=None,
    ts_group="TS", ct_group="CT",
    ts_exclude_modules=None,
    ct_exclude_modules=None,
    hist_as_proportion = False
):
    assert build_matrix_fn is not None, "אנא ספק את build_kegg_logp_matrix דרך הפרמטר build_matrix_fn"
    if show_region_legend is None:
        show_region_legend = bool(add_hist and hist_stacked_regions)
    regions = [str(r).strip() for r in hist_regions_cols]
    if hist_region_colors is None:
        hist_region_colors = {"AC":"#88CCEE","MF":"#CC6677","PCG":"#DDCC77"}
    default_pool = list(plt.cm.Set2.colors)
    region_color_map = {}
    for i, r in enumerate(regions):
        region_color_map[r] = hist_region_colors.get(r, default_pool[i % len(default_pool)])

    logp_ts, meas_ts, stars_ts = build_matrix_fn(
        kegg_long_csv, modules_details=modules_details_tsv, clean_TS=True,
        cap=cap, selection=selection, K_GLOBAL=K_GLOBAL,
        N_PER_MODULE=N_PER_MODULE, P_THRESH=P_THRESH,
        star_from=star_from, group=ts_group
    )
    logp_ct, meas_ct, stars_ct = build_matrix_fn(
        kegg_long_csv, modules_details=modules_details_tsv, clean_TS=True,
        cap=cap, selection=selection, K_GLOBAL=K_GLOBAL,
        N_PER_MODULE=N_PER_MODULE, P_THRESH=P_THRESH,
        star_from=star_from, group=ct_group
    )

    def _drop_cols(df, drop_set):
        if not drop_set: return df
        keep = []
        for c in df.columns:
            try:
                ci = int(c)
            except Exception:
                try:
                    ci = int(str(c).strip())
                except Exception:
                    ci = None
            if (ci is None) or (ci not in drop_set):
                keep.append(c)
        return df[keep] if len(keep) else df.iloc[:, :0]

    ts_drop = _normalize_mod_set(ts_exclude_modules)
    ct_drop = _normalize_mod_set(ct_exclude_modules)
    if len(ts_drop) > 0:
        logp_ts   = _drop_cols(logp_ts, ts_drop)
        if stars_ts is not None: stars_ts = stars_ts.reindex(index=logp_ts.index, columns=logp_ts.columns)
    if len(ct_drop) > 0:
        logp_ct   = _drop_cols(logp_ct, ct_drop)
        if stars_ct is not None: stars_ct = stars_ct.reindex(index=logp_ct.index, columns=logp_ct.columns)

    row_meta_ts = _extract_row_meta_from_kegg(kegg_long_csv, logp_ts.index.tolist())
    row_meta_ct = _extract_row_meta_from_kegg(kegg_long_csv, logp_ct.index.tolist())
    if palette_cat is None or palette_sub is None:
        palette_cat_all, palette_sub_all = build_annotation_palettes_from_kegg(kegg_long_csv)
        if palette_cat is None: palette_cat = palette_cat_all
        if palette_sub is None: palette_sub = palette_sub_all
    cat_order, sub_order = _appearance_orders_from_kegg(kegg_long_csv)

    ts_row_order, ts_col_order, ts_seps = _cluster_orders(
        logp_ts, row_meta_ts, cluster_cols, col_metric, col_method,
        cluster_rows, row_metric, row_method, row_grouping=row_grouping,
        cat_order=cat_order, sub_order=sub_order
    )
    ct_row_order, ct_col_order, ct_seps = _cluster_orders(
        logp_ct, row_meta_ct, cluster_cols, col_metric, col_method,
        cluster_rows, row_metric, row_method, row_grouping=row_grouping,
        cat_order=cat_order, sub_order=sub_order
    )

    logp_ts = logp_ts.loc[ts_row_order, ts_col_order]
    logp_ct = logp_ct.loc[ct_row_order, ct_col_order]
    if stars_ts is not None: stars_ts = stars_ts.reindex(index=ts_row_order, columns=ts_col_order)
    if stars_ct is not None: stars_ct = stars_ct.reindex(index=ct_row_order, columns=ct_col_order)
    row_meta_ts = row_meta_ts.loc[ts_row_order]
    row_meta_ct = row_meta_ct.loc[ct_row_order]

    def _arr_from_palette(labels, pal, default=(0.92,0.92,0.92)):
        def lk(v):
            if pd.isna(v): return default
            vv = str(v).strip()
            return pal.get(vv, default)
        return np.array([to_rgb(lk(v)) for v in labels])

    ts_cat_colors = _arr_from_palette(row_meta_ts["category"], palette_cat)
    ts_sub_colors = _arr_from_palette(row_meta_ts["subcategory"], palette_sub)
    ct_cat_colors = _arr_from_palette(row_meta_ct["category"], palette_cat)
    ct_sub_colors = _arr_from_palette(row_meta_ct["subcategory"], palette_sub)

    def _canon(x): return str(x).strip().upper()
    regions = [_canon(r) for r in hist_regions_cols]
    region_df = _read_region_counts(modules_details_tsv, regions).rename(columns=_canon)




    mod_sizes = _read_module_sizes(modules_details_tsv) if add_hist and not hist_stacked_regions else None

    strip = strip_width; gap = 0.05; legend_w = 0.35
    hr_hist = (0.18 if add_hist else 0.01); hr_hm, hr_cbar = 1.00, 0.02

    fig = plt.figure(figsize=(22, 12), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=3, ncols=8,
        width_ratios=[strip, strip, 1.0, gap, strip, strip, 1.0, legend_w],
        height_ratios=[hr_hist, hr_hm, hr_cbar],
        wspace=0.02, hspace=0.02
    )

    ax_ts_hist = fig.add_subplot(gs[0,2]) if add_hist else fig.add_subplot(gs[0,2], frame_on=False)
    ax_ts_sub  = fig.add_subplot(gs[1,0]); ax_ts_cat = fig.add_subplot(gs[1,1]); ax_ts_hm = fig.add_subplot(gs[1,2])
    ax_ts_cbar = fig.add_subplot(gs[2,2])

    ax_ct_hist = fig.add_subplot(gs[0,6]) if add_hist else fig.add_subplot(gs[0,6], frame_on=False)
    ax_ct_sub  = fig.add_subplot(gs[1,4]); ax_ct_cat = fig.add_subplot(gs[1,5]); ax_ct_hm = fig.add_subplot(gs[1,6])
    ax_ct_cbar = fig.add_subplot(gs[2,6])

    ax_leg = fig.add_subplot(gs[:,7]); ax_leg.axis("off")
    ax_ts_hist.sharex(ax_ts_hm)
    ax_ct_hist.sharex(ax_ct_hm)

    ax_ts_hist.set_xlim(ax_ts_hm.get_xlim())
    ax_ct_hist.set_xlim(ax_ct_hm.get_xlim())
    data_ts = np.nan_to_num(logp_ts.values, nan=0.0)
    im_ts = ax_ts_hm.imshow(data_ts, aspect="auto", vmin=0.0, vmax=cap, cmap="Blues")
    cb_ts = fig.colorbar(im_ts, cax=ax_ts_cbar, orientation="horizontal")
    cb_ts.ax.tick_params(labelsize=8); cb_ts.set_label(f"−log10(p) (0..{cap})", labelpad=5)
    ax_ts_hm.set_xticks(range(logp_ts.shape[1])); ax_ts_hm.set_xticklabels(logp_ts.columns.astype(str), rotation=90)
    ax_ts_hm.set_yticks([])


    if row_grouping is not None and len(ts_seps) > 0:
        for y in ts_seps:
            ax_ts_hm.hlines(y, -0.5, logp_ts.shape[1]-0.5, colors="k", linewidth=0.4, alpha=0.4)

    ax_ts_sub.imshow(ts_sub_colors.reshape(-1,1,3), aspect="auto")
    ax_ts_cat.imshow(ts_cat_colors.reshape(-1,1,3), aspect="auto")
    for ax in (ax_ts_sub, ax_ts_cat):
        ax.set_ylim(ax_ts_hm.get_ylim()); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
    # --- כוכביות TS ---
    if stars_ts is not None and stars_ts.size:
        S = stars_ts.reindex(index=logp_ts.index, columns=logp_ts.columns)
        for i, row in enumerate(S.values):
            for j, sym in enumerate(row):
                if isinstance(sym, str) and sym:
                    ax_ts_hm.text(
                        j, i, sym, ha="center", va="center",
                        color=star_color, fontsize=star_fontsize, zorder=10,
                        path_effects=[pe.Stroke(linewidth=1.4, foreground="white"), pe.Normal()] if star_outline else None
                    )

    if add_hist:
        x = np.arange(len(logp_ts.columns))
        ts_mod_ids = pd.Index(logp_ts.columns).astype(int)

        if hist_stacked_regions and region_df is not None:
            cols_int = pd.Index(logp_ts.columns).astype(int)
            rsel_ts = region_df.reindex(index=ts_mod_ids, columns=regions).fillna(0)
            assert np.array_equal(rsel_ts.index.values, ts_mod_ids.values)
            if hist_as_proportion:
                rsel_ts = rsel_ts.div(rsel_ts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
            x_ts = np.arange(len(ts_mod_ids))
            bottom = np.zeros(len(x_ts))
            for rc in regions:   
                vals = rsel_ts[rc].to_numpy() if rc in rsel_ts.columns else np.zeros(len(x_ts))
                ax_ts_hist.bar(x_ts, vals, width=bar_width, color=region_color_map[rc], bottom=bottom, label=rc)
                bottom += vals

            ax_ts_hist.set_xlim(-0.5, len(x_ts)-0.5); ax_ts_hist.set_xticks([])
            if hist_as_proportion:
                ax_ts_hist.set_ylim(0, 1)
                ax_ts_hist.set_yticks([0, 0.5, 1.0])
                ax_ts_hist.set_yticklabels(["0%", "50%", "100%"])
                ax_ts_hist.set_ylabel("Region %", fontsize=9)
            else:
                ax_ts_hist.set_ylabel("Genes", fontsize=9)

        else:
            sizes = (mod_sizes.reindex(pd.Index(logp_ts.columns).astype(int)).fillna(0).values
                    if mod_sizes is not None else np.zeros(len(x)))
            ax_ts_hist.bar(x, sizes, width=bar_width, color=bar_color)
            ax_ts_hist.set_xlim(-0.5, len(x)-0.5); ax_ts_hist.set_xticks([])
            ax_ts_hist.set_ylabel("Genes", fontsize=9)
    ax_ts_hist.set_title("TS KEGG enrichment")

    data_ct = np.nan_to_num(logp_ct.values, nan=0.0)
    im_ct = ax_ct_hm.imshow(data_ct, aspect="auto", vmin=0.0, vmax=cap, cmap="Blues")
    cb_ct = fig.colorbar(im_ct, cax=ax_ct_cbar, orientation="horizontal")
    cb_ct.ax.tick_params(labelsize=8); cb_ct.set_label(f"−log10(p) (0..{cap})", labelpad=5)
    ax_ct_hm.set_xticks(range(logp_ct.shape[1])); ax_ct_hm.set_xticklabels(logp_ct.columns.astype(str), rotation=90)
    ax_ct_hm.set_yticks([])
    

    if row_grouping is not None and len(ct_seps) > 0:
        for y in ct_seps:
            ax_ct_hm.hlines(y, -0.5, logp_ct.shape[1]-0.5, colors="k", linewidth=0.4, alpha=0.4)

    ax_ct_sub.imshow(ct_sub_colors.reshape(-1,1,3), aspect="auto")
    ax_ct_cat.imshow(ct_cat_colors.reshape(-1,1,3), aspect="auto")
    for ax in (ax_ct_sub, ax_ct_cat):
        ax.set_ylim(ax_ct_hm.get_ylim()); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
    # --- כוכביות CT ---
    if stars_ct is not None and stars_ct.size:
        S = stars_ct.reindex(index=logp_ct.index, columns=logp_ct.columns)
        for i, row in enumerate(S.values):
            for j, sym in enumerate(row):
                if isinstance(sym, str) and sym:
                    ax_ct_hm.text(
                        j, i, sym, ha="center", va="center",
                        color=star_color, fontsize=star_fontsize, zorder=10,
                        path_effects=[pe.Stroke(linewidth=1.4, foreground="white"), pe.Normal()] if star_outline else None
                    )

    if add_hist:
        ct_mod_ids = pd.Index(logp_ct.columns).astype(int)
        x = np.arange(len(logp_ct.columns))
        x_ct = np.arange(len(ct_mod_ids))
        if hist_stacked_regions and region_df is not None:
            cols_int = pd.Index(logp_ct.columns).astype(int)
            rsel_ct = region_df.reindex(index=ct_mod_ids, columns=regions).fillna(0)
            assert np.array_equal(rsel_ct.index.values, cols_int.values), \
                f"Misalignment CT: expected {list(cols_int[:5])}, got {list(rsel_ct.index[:5])}"
            if hist_as_proportion:
                totals = rsel_ct.sum(axis=1).replace(0, np.nan)
                rsel_ct = rsel_ct.div(rsel_ct.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)

            bottom = np.zeros(len(x_ct))
            for rc in regions:
                vals = rsel_ct[rc].to_numpy()
                ax_ct_hist.bar(x_ct, vals, width=bar_width, color=region_color_map[rc], bottom=bottom, label=rc)
                bottom += vals

            ax_ct_hist.set_xlim(-0.5, len(x_ct)-0.5); ax_ct_hist.set_xticks([])
            if hist_as_proportion:
                ax_ct_hist.set_ylim(0, 1)
                ax_ct_hist.set_yticks([0, 0.5, 1.0])
                ax_ct_hist.set_yticklabels(["0%", "50%", "100%"])
                ax_ct_hist.set_ylabel("Region %", fontsize=9)
                for j, mid in enumerate(ct_mod_ids):
                    ax_ct_hist.text(j, ax_ct_hist.get_ylim()[1]*0.98, str(mid),
                                    ha="center", va="top", fontsize=6, rotation=90, alpha=0.6)
            else:
                ax_ct_hist.set_ylabel("Genes", fontsize=9)

        else:
            sizes = (mod_sizes.reindex(pd.Index(logp_ct.columns).astype(int)).fillna(0).values
                    if mod_sizes is not None else np.zeros(len(x)))
            ax_ct_hist.bar(x, sizes, width=bar_width, color=bar_color)
            ax_ct_hist.set_xlim(-0.5, len(x)-0.5); ax_ct_hist.set_xticks([])
            ax_ct_hist.set_ylabel("Genes", fontsize=9)
    ax_ct_hist.set_title("CT KEGG enrichment", pad=2, fontsize=11)
    handles_cat = [Patch(facecolor=palette_cat[k], edgecolor='none', label=str(k)) for k in palette_cat]
    handles_sub = [Patch(facecolor=palette_sub[k], edgecolor='none', label=str(k)) for k in palette_sub]
    
    def place_leg(ax, handles, title, y_top, fs):
        if not handles: return y_top
        est_h = 0.035 * (len(handles) + 2)
        leg = ax.legend(handles=handles, title=title, loc="upper left",
                        bbox_to_anchor=(0.0, y_top), fontsize=fs, title_fontsize=fs+1, frameon=False)
        ax.add_artist(leg)
        return max(0.02, y_top - est_h)

    y = 0.98
    y = place_leg(ax_leg, handles_cat,     "Category",     y, legend_fontsize)
    y = place_leg(ax_leg, handles_sub,     "Subcategory",  y, legend_fontsize)
    if show_region_legend:
        handles_regions = [Patch(facecolor=region_color_map[r], edgecolor='none', label=str(r)) for r in regions]
        _ = place_leg(ax_leg, handles_regions, "Regions (hist)", y, legend_fontsize)


    plt.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.06, wspace=0.02, hspace=0.02)
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    return fig, ts_row_order, ts_col_order, ct_row_order, ct_col_order

def build_kegg_logp_matrix(
    kegg_long_csv,           
    modules_details=None,   
    clean_TS=True,
    cap=6.0,                 
    selection="global_top", 
    K_GLOBAL=60,             
    N_PER_MODULE=4,             
    P_THRESH=0.05,              
    star_from="p",
    group = 'TS',
    alpha_levels=((0.0005,"***"), (0.005,"**"), (0.05,"*"))
):
    df = pd.read_csv(kegg_long_csv, engine="python", on_bad_lines="skip")
    assert "Cluster" in df.columns and "Description" in df.columns, "Missing Cluster/Description"

    df["Cluster"] = pd.to_numeric(df["Cluster"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["Cluster"])
    df["Cluster"] = df["Cluster"].astype(int)
    df["Description"] = df["Description"].astype(str).str.strip()
    if "ID" in df.columns:
        df["ID"] = df["ID"].astype(str).str.strip()

    df["p"]      = pd.to_numeric(df.get("pvalue"),   errors="coerce")
    df["padj"]   = pd.to_numeric(df.get("p.adjust"), errors="coerce")
    df["qvalue"] = pd.to_numeric(df.get("qvalue"),   errors="coerce")

    if clean_TS and modules_details is not None:
        md = pd.read_csv(modules_details, sep="\t")
        ct_ids = md.loc[md["Cluster Type"].eq(group), "Cluster ID"].astype(int).tolist()
        df = df[df["Cluster"].isin(ct_ids)]
    

    agg = (df.groupby(["Description","Cluster"], as_index=False)
             .agg(p=("p","min"), padj=("padj","min"), qvalue=("qvalue","min")))
    pmat = agg.pivot(index="Description", columns="Cluster", values="p").astype(float)

    if selection == "global_top":
        chosen_terms = pmat.min(axis=1).nsmallest(K_GLOBAL).index.tolist()
    elif selection == "per_module":
        chosen = set()
        for c in pmat.columns:
            col = pmat[c]
            sig = col[col < P_THRESH].nsmallest(N_PER_MODULE)
            if sig.empty: sig = col.nsmallest(N_PER_MODULE)
            chosen.update(sig.index.tolist())
        chosen_terms = list(chosen)
        if len(chosen_terms) > K_GLOBAL:
            chosen_terms = pmat.loc[chosen_terms].min(axis=1).nsmallest(K_GLOBAL).index.tolist()
    else:
        raise ValueError("selection must be 'global_top' או 'per_module'.")

    psub = pmat.loc[chosen_terms] 

    EPS = np.nextafter(0, 1)
    logp_df = -np.log10(np.clip(psub, EPS, 1.0))
    logp_df = logp_df.clip(lower=0, upper=cap).fillna(0.0)

    if "ID" in df.columns:
        idmap = (df[["Description","ID"]].dropna().drop_duplicates()
                   .groupby("Description")["ID"].first())
        new_index = [f"{idmap.get(desc, '')} — {desc}".strip(" —") for desc in logp_df.index]
        logp_df.index = new_index

    if star_from == "padj":
        measure_mat = agg.pivot(index="Description", columns="Cluster", values="padj").astype(float)
    elif star_from in {"qvalue","q"}:
        measure_mat = agg.pivot(index="Description", columns="Cluster", values="qvalue").astype(float)
    else:
        measure_mat = pmat  

    measure_df = measure_mat.loc[psub.index, psub.columns]
    if "ID" in df.columns:
        measure_df.index = logp_df.index  

    def _stars(p):
        if not np.isfinite(p): return ""
        for thr, sym in alpha_levels:
            if p <= thr: return sym
        return ""
    stars_df = measure_df.applymap(_stars)

    max_val = logp_df.max(axis=1)
    peak_mod = logp_df.idxmax(axis=1)
    logp_df["__peak_mod__"] = peak_mod
    logp_df["__max__"] = max_val
    logp_df = logp_df.sort_values(by=["__peak_mod__", "__max__"], ascending=[True, False]) \
                     .drop(columns=["__peak_mod__", "__max__"])

    measure_df = measure_df.loc[logp_df.index, logp_df.columns]
    stars_df   = stars_df.loc[logp_df.index, logp_df.columns]

    return logp_df, measure_df, stars_df 


def plot_log_heatmap_with_stars(logp_df, stars_df=None, cap=6.0, title="KEGG −log10(p)"):
    fig, ax = plt.subplots(figsize=(18, 10))
    im = ax.imshow(logp_df.values, aspect="auto", vmin=0.0, vmax=cap)
    fig.colorbar(im, ax=ax, label=f"−log10(p) (0..{cap})")
    ax.set_yticks(range(logp_df.shape[0])); ax.set_yticklabels(logp_df.index, fontsize=8)
    ax.set_xticks(range(logp_df.shape[1])); ax.set_xticklabels(logp_df.columns.astype(str), rotation=90)
    ax.set_xlabel("Module"); ax.set_title(title)

    if stars_df is not None:
        stars_df = stars_df.reindex(index=logp_df.index, columns=logp_df.columns)
        for i, term in enumerate(logp_df.index):
            for j, mod in enumerate(logp_df.columns):
                s = stars_df.loc[term, mod]
                if isinstance(s, str) and s:
                    ax.text(j, i, s, ha="center", va="center", fontsize=9, color="black")

    fig.tight_layout()
    return fig



def _normalize_mod_set(lst):
    if lst is None: return set()
    norm = set()
    for x in lst:
        try:
            norm.add(int(x))
        except Exception:
            try:
                norm.add(int(str(x).strip()))
            except Exception:
                pass
    return norm




def _normalize_mod_set(lst):
    if lst is None: return set()
    norm = set()
    for x in lst:
        try:
            norm.add(int(x))
        except Exception:
            try:
                norm.add(int(str(x).strip()))
            except Exception:
                pass
    return norm


def _read_region_counts(modules_details_tsv, region_cols=("AC","MF","PCG")):
    md  = pd.read_csv(modules_details_tsv, sep="\t")
    cid = pd.to_numeric(md["Cluster ID"], errors="coerce").astype("Int64")
    md  = md.loc[cid.notna()].copy()
    cid = cid.astype(int)

    data = {
        c: (pd.to_numeric(md[c], errors="coerce").fillna(0).astype(int).to_numpy()
            if c in md.columns else np.zeros(len(md), dtype=int))
        for c in region_cols
    }
    df = pd.DataFrame(data)         
    df.index = cid.to_numpy()         
    df = df.groupby(level=0, sort=True).sum()
    df.index = df.index.astype(int)
    return df

def _normalize_mod_set(lst):
    if lst is None: return set()
    norm = set()
    for x in lst:
        try:
            norm.add(int(x))
        except Exception:
            try:
                norm.add(int(str(x).strip()))
            except Exception:
                pass
    return norm

def _appearance_orders_from_kegg(kegg_long_csv):
    df = pd.read_csv(kegg_long_csv, engine="python", on_bad_lines="skip")
    cols = _std_cols(df)
    def _clean_series(s):
        if s is None or s not in df.columns:
            return []
        ss = df[s].astype(str).str.strip()
        uniq = [x for x in pd.unique(ss) if x and x.lower() not in ("nan","none")]
        return uniq
    cat_order = _clean_series(cols["category"])
    sub_order = _clean_series(cols["subcategory"])
    return cat_order, sub_order
def _cluster_orders(logp_df, row_meta,
                    cluster_cols=True, col_metric="euclidean", col_method="average",
                    cluster_rows=True, row_metric="euclidean", row_method="average",
                    row_grouping=("category","subcategory"),
                    cat_order=None, sub_order=None):
    if cluster_cols:
        gcol = sns.clustermap(
            np.nan_to_num(logp_df.values, nan=0.0),
            row_cluster=False, col_cluster=True,
            metric=col_metric, method=col_method,
            cmap="Blues", cbar=False, xticklabels=False, yticklabels=False
        )
        col_order = [logp_df.columns[i] for i in gcol.dendrogram_col.reordered_ind]
        plt.close(gcol.fig)
    else:
        col_order = list(logp_df.columns)

    rm = row_meta.copy()
    def _norm(x): return np.nan if pd.isna(x) else str(x).strip()
    rm["category"]    = rm["category"].map(_norm)
    rm["subcategory"] = rm["subcategory"].map(_norm)

    sep_positions = []

    if not cluster_rows:
        row_order = list(logp_df.index)
    else:
        if row_grouping is None:
            grow = sns.clustermap(
                np.nan_to_num(logp_df.values, nan=0.0),
                row_cluster=True, col_cluster=False,
                metric=row_metric, method=row_method,
                cmap="Blues", cbar=False, xticklabels=False, yticklabels=False
            )
            row_order = [logp_df.index[i] for i in grow.dendrogram_row.reordered_ind]
            plt.close(grow.fig)
        else:
            def _cats_for(col_name, global_order):
                seen_local = [x for x in pd.unique(rm[col_name].dropna())]
                if global_order:
                    return list(global_order) + [x for x in seen_local if x not in global_order]
                return seen_local

            cat_cats = _cats_for("category",    cat_order)
            sub_cats = _cats_for("subcategory", sub_order)

            keys = list(row_grouping)
            for k in keys:
                if k == "category":
                    cats = cat_cats
                elif k == "subcategory":
                    cats = sub_cats
                else:
                    cats = [x for x in pd.unique(rm[k].dropna())]
                rm[f"__ord_{k}__"] = pd.Categorical(rm[k], categories=cats, ordered=True)

            sorted_rm = rm.sort_values(by=[f"__ord_{k}__" for k in keys],
                                       na_position="last", kind="mergesort")

            row_order = []
            last_cat = None; cum_len = 0
            for grp_keys, sub_rm in sorted_rm.groupby(keys, sort=False):
                idxs = sub_rm.index
                sub  = logp_df.loc[idxs]
                if len(sub) > 1:
                    g = sns.clustermap(
                        np.nan_to_num(sub.values, nan=0.0),
                        row_cluster=True, col_cluster=False,
                        metric=row_metric, method=row_method,
                        cmap="Blues", cbar=False, xticklabels=False, yticklabels=False
                    )
                    reord = [sub.index[i] for i in g.dendrogram_row.reordered_ind]
                    plt.close(g.fig)
                else:
                    reord = list(idxs)
                curr_cat = sub_rm["category"].iloc[0]
                if last_cat is None:
                    last_cat = curr_cat
                elif curr_cat != last_cat:
                    sep_positions.append(cum_len - 0.5)
                    last_cat = curr_cat
                row_order.extend(reord)
                cum_len += len(reord)

    return row_order, col_order, sep_positions



if __name__ == "__main__":

    KEGG_CSV = r"/Users/edeneldar/Downloads/kegg_rosmap_full_V4.csv"

    MODULES = "/Users/edeneldar/Downloads/tryRosmap_Cluster_details3.tsv"

    palette_cat, palette_sub = build_annotation_palettes_from_kegg(KEGG_CSV)

    fig, ts_rows, ts_cols, ct_rows, ct_cols = export_kegg_ts_ct_enrichment_pdf(
        kegg_long_csv=KEGG_CSV,
        modules_details_tsv=MODULES,
        pdf_path="kegg_TS_CT_enrichment_between_V4.pdf",
        selection="global_top", K_GLOBAL=60, cap=4.0,
        cluster_cols=True, cluster_rows=True, row_grouping=(None),
        add_hist=True,
        hist_stacked_regions=True,
        hist_regions_cols=(["AC","MFBA9BA46","PCGBA23"]),
        hist_region_colors={"AC":"#88CCEE","MFBA9BA46":"#CC6677","PCGBA23":"#DDCC77"},
        palette_cat=palette_cat, palette_sub=palette_sub,
        build_matrix_fn=build_kegg_logp_matrix,
        ts_exclude_modules=[],           
        ct_exclude_modules=[],
        star_color="crimson", star_fontsize=6, star_outline=False, hist_as_proportion=True, show_row_labels=True
    )
