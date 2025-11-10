#!/usr/bin/env python3
"""
Fetch a small scRNA-seq dataset using Scanpy's dataset loader.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import scanpy as sc


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-file", required=True)
    parser.add_argument("--summary-json", required=True)
    args = parser.parse_args(list(argv) if argv is not None else None)

    adata = sc.datasets.pbmc3k_processed()
    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)

    summary = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "obs_keys": list(adata.obs.columns),
        "var_keys": list(adata.var.columns),
    }
    Path(args.summary_json).write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

