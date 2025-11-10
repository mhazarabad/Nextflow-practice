#!/usr/bin/env python3
"""
Summarise scRNA-seq expression of a gene set across clusters.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata", required=True, help="Path to .h5ad file.")
    parser.add_argument("--gene-table", required=True, help="Gene list TSV with gene_name column.")
    parser.add_argument("--cluster-key", default="louvain")
    parser.add_argument("--out-signature", required=True)
    parser.add_argument("--summary-json", required=True)
    args = parser.parse_args(list(argv) if argv is not None else None)

    adata = ad.read_h5ad(args.adata)
    if args.cluster_key not in adata.obs.columns:
        raise RuntimeError(f"Cluster key '{args.cluster_key}' not present in obs.")

    gene_table = pd.read_csv(args.gene_table, sep="\t")
    if "gene_name" not in gene_table.columns:
        raise RuntimeError("Gene table must include 'gene_name' column.")

    target_genes = sorted(set(gene_table["gene_name"]).intersection(set(adata.var_names)))
    if not target_genes:
        raise RuntimeError("No overlap between gene list and scRNA dataset.")

    # ensure dense matrix for simple operations
    matrix = adata[:, target_genes].X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix)

    df = pd.DataFrame(matrix, columns=target_genes, index=adata.obs[args.cluster_key])
    grouped = df.groupby(level=0)
    mean_expr = grouped.mean()
    median_expr = grouped.median()

    signature = mean_expr.stack().reset_index()
    signature.columns = ["cluster", "gene_name", "mean_expression"]
    signature["median_expression"] = median_expr.stack().values
    signature.to_csv(args.out_signature, sep="\t", index=False)

    # summarise top cluster per gene
    top_clusters = (
        signature.groupby("gene_name")
        .apply(lambda x: x.sort_values("mean_expression", ascending=False).iloc[0][["cluster", "mean_expression"]])
        .reset_index()
        .rename(columns={"cluster": "top_cluster", "mean_expression": "top_cluster_mean"})
    )

    summary = {
        "n_clusters": int(signature["cluster"].nunique()),
        "n_genes_profiled": int(len(target_genes)),
        "top_cluster_assignments": top_clusters.head(10).to_dict(orient="records"),
    }
    Path(args.summary_json).write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

