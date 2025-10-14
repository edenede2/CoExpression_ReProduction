#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GTEx Packager: collect young/old matrices per tissue and build full (vertical concat).

Expected input layout:
Tissues/
  <TISSUE_A>/Raw/young_matrix.csv
  <TISSUE_A>/Raw/old_matrix.csv
  <TISSUE_B>/Raw/young_matrix.csv
  <TISSUE_B>/Raw/old_matrix.csv
  ...

Outputs (under --out root, default: <Tissues>/GTEX_outputs):
  GTEX_outputs/young/GTEX_young_<tissueStrippedName>.csv
  GTEX_outputs/old/  GTEX_old_<tissueStrippedName>.csv
  GTEX_outputs/full/ GTEX_full_<tissueStrippedName>.csv   (young stacked over old)

Notes
- tissueStrippedName is derived from the tissue folder name:
  * --strategy camel (default):  "Skin - Sun Exposed (Lower leg)" -> "SkinSunExposedLowerLeg"
  * --strategy alnum:            keep only [A-Za-z0-9], joined     -> "SkinSunExposedLowerleg"
  * --strategy snake:            lower + '_'                       -> "skin_sun_exposed_lower_leg"
  * --strategy kebab:            lower + '-'                       -> "skin-sun-exposed-lower-leg"
- If one of young/old is missing, the script still writes the available single file and also uses it for "full".
- Column mismatch is handled by outer-union of columns (missing columns filled with nulls).
- Requires Python 3.9+. Prefers polars; falls back to pandas if polars isn't installed.

Usage
------
python gtex_packager.py /path/to/Tissues
# with options
python gtex_packager.py /path/to/Tissues --out /path/to/out --strategy alnum --dry-run
"""

import argparse
import sys
import re
from pathlib import Path
from typing import List, Tuple, Optional

def sanitize_tissue(name: str, strategy: str = "camel") -> str:
    tokens = [t for t in re.split(r"[^A-Za-z0-9]+", name) if t]
    if strategy == "camel":
        return "".join(t[0].upper() + t[1:].lower() if t else "" for t in tokens)
    elif strategy == "snake":
        return "_".join(t.lower() for t in tokens)
    elif strategy == "kebab":
        return "-".join(t.lower() for t in tokens)
    elif strategy == "alnum":
        return "".join(ch for ch in name if ch.isalnum())
    else:
        return "".join(tokens)

def ensure_dir(p: Path, dry_run: bool = False) -> None:
    if dry_run:
        print(f"[dry-run] mkdir -p {p}")
        return
    p.mkdir(parents=True, exist_ok=True)

def find_tissue_dirs(tissues_root: Path) -> List[Path]:
    candidates = sorted([p for p in tissues_root.iterdir() if p.is_dir() and not p.name.startswith(".")])
    tissue_dirs = [p for p in candidates if (p / "Raw").is_dir()]
    return tissue_dirs

def load_csv(path: Path):
    try:
        import polars as pl  # type: ignore
        df = pl.read_csv(str(path))
        return ("polars", df)
    except Exception as e_pl:
        try:
            import pandas as pd  # type: ignore
            df = pd.read_csv(path)
            return ("pandas", df)
        except Exception as e_pd:
            raise RuntimeError(f"Failed to read CSV '{path}':\n - polars error: {e_pl}\n - pandas error: {e_pd}")

def save_csv(path: Path, backend: str, df, dry_run: bool = False) -> None:
    if dry_run:
        print(f"[dry-run] write {path}")
        return
    if backend == "polars":
        df.write_csv(str(path))
    else:
        df.to_csv(path, index=False)

def union_columns(backend: str, a, b):
    if backend == "polars":
        import polars as pl  # type: ignore
        cols_a = set(a.columns)
        cols_b = set(b.columns)
        all_cols = list(cols_a | cols_b)
        def add_missing(df, all_cols: List[str]):
            exprs = []
            for c in all_cols:
                if c in df.columns:
                    exprs.append(pl.col(c))
                else:
                    exprs.append(pl.lit(None).alias(c))
            return df.select(exprs)
        a2 = add_missing(a, all_cols)
        b2 = add_missing(b, all_cols)
        return a2, b2, all_cols
    else:
        import pandas as pd  # type: ignore
        all_cols = list(set(a.columns) | set(b.columns))
        a2 = a.reindex(columns=all_cols)
        b2 = b.reindex(columns=all_cols)
        return a2, b2, all_cols

def vconcat(backend: str, frames: List):
    if backend == "polars":
        import polars as pl  # type: ignore
        return pl.concat(frames, how="vertical", rechunk=True)
    else:
        import pandas as pd  # type: ignore
        return pd.concat(frames, axis=0, ignore_index=True)

def process_tissue(tissue_dir: Path, out_root: Path, strategy: str, dry_run: bool, overwrite: bool) -> Tuple[bool, str]:
    raw_dir = tissue_dir / "Raw"
    y_path = raw_dir / "young_matrix.csv"
    o_path = raw_dir / "old_matrix.csv"

    if not y_path.exists() and not o_path.exists():
        return False, f"SKIP: no young/old CSV found under {raw_dir}"

    tissue_name = tissue_dir.name
    stripped = sanitize_tissue(tissue_name, strategy=strategy)

    out_y = out_root / "young" / f"GTEX_young_{stripped}.csv"
    out_o = out_root / "old"  / f"GTEX_old_{stripped}.csv"
    out_f = out_root / "full" / f"GTEX_full_{stripped}.csv"

    ensure_dir(out_y.parent, dry_run=dry_run)
    ensure_dir(out_o.parent, dry_run=dry_run)
    ensure_dir(out_f.parent, dry_run=dry_run)

    backend = None
    y_df = o_df = None

    if y_path.exists():
        backend, y_df = load_csv(y_path)
    if o_path.exists():
        if backend is None:
            backend, o_df = load_csv(o_path)
        else:
            try:
                if backend == "polars":
                    import polars as pl  # type: ignore
                    o_df = pl.read_csv(str(o_path))
                else:
                    import pandas as pd  # type: ignore
                    o_df = pd.read_csv(o_path)
            except Exception:
                _, o_df = load_csv(o_path)

    if y_df is not None:
        if not overwrite and out_y.exists():
            print(f"SKIP (exists): {out_y}")
        else:
            save_csv(out_y, backend, y_df, dry_run=dry_run)

    if o_df is not None:
        if not overwrite and out_o.exists():
            print(f"SKIP (exists): {out_o}")
        else:
            save_csv(out_o, backend, o_df, dry_run=dry_run)

    if y_df is not None and o_df is not None:
        a2, b2, _ = union_columns(backend, y_df, o_df)
        full_df = vconcat(backend, [a2, b2])
    elif y_df is not None:
        full_df = y_df
    elif o_df is not None:
        full_df = o_df
    else:
        return False, f"SKIP: neither young nor old were loadable for {tissue_dir}"

    if not overwrite and out_f.exists():
        print(f"SKIP (exists): {out_f}")
    else:
        save_csv(out_f, backend, full_df, dry_run=dry_run)

    return True, f"DONE: {tissue_name} -> {stripped}"

def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Pack GTEx per-tissue young/old matrices into three folders and build full CSVs.")
    p.add_argument("tissues_dir", type=Path, help="Path to Tissues/ root directory")
    p.add_argument("--out", type=Path, default=None, help="Output root directory (default: <Tissues>/GTEX_outputs)")
    p.add_argument("--strategy", choices=["camel", "alnum", "snake", "kebab"], default="camel", help="How to strip tissue names (default: camel)")
    p.add_argument("--dry-run", action="store_true", help="Print actions without writing files")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    args = p.parse_args(argv)

    tissues_root = args.tissues_dir
    if not tissues_root.exists() or not tissues_root.is_dir():
        print(f"ERROR: '{tissues_root}' is not a directory.", file=sys.stderr)
        return 2

    out_root = args.out if args.out is not None else (tissues_root / "GTEX_outputs")
    if not args.dry_run:
        ensure_dir(out_root, dry_run=False)

    tissue_dirs = find_tissue_dirs(tissues_root)
    if not tissue_dirs:
        print("No tissue directories with a 'Raw' subfolder were found. Nothing to do.")
        return 0

    print(f"Discovered {len(tissue_dirs)} tissue(s). Writing to: {out_root}")
    print(f"Name strategy: {args.strategy}")

    ok = 0
    skipped = 0
    for td in tissue_dirs:
        success, msg = process_tissue(td, out_root, args.strategy, args.dry_run, args.overwrite)
        print(msg)
        if success:
            ok += 1
        else:
            skipped += 1

    print(f"\nSummary: processed={ok}, skipped={skipped}, total={len(tissue_dirs)}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
