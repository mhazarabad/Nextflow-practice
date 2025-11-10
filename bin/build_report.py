#!/usr/bin/env python3
"""
Integrate bulk RNA-seq DEGs, ChIP targets, and scRNA signatures.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--deg", required=True)
    parser.add_argument("--chip-genes", required=True)
    parser.add_argument("--scrna-signature", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--top-n", type=int, default=25)
    args = parser.parse_args(list(argv) if argv is not None else None)

    deg = pd.read_csv(args.deg, sep="\t")
    chip = pd.read_csv(args.chip_genes, sep="\t")
    scrna = pd.read_csv(args.scrna_signature, sep="\t")

    merged = (
        deg.merge(chip, how="inner", on=["gene_id", "gene_name"])
        .merge(
            scrna,
            how="left",
            on="gene_name",
            suffixes=("", "_scrna"),
        )
        .sort_values(["padj", "n_peaks", "mean_expression"], ascending=[True, False, False])
    )

    merged.to_csv(args.out_summary, sep="\t", index=False)

    figure_dir = Path(args.out_summary).parent
    heatmap_path = figure_dir / "multiomic_heatmap.png"
    bar_path = figure_dir / "chip_target_bar.png"

    _plot_heatmap(merged, heatmap_path)
    _plot_bar(merged, bar_path)

    # build compact summary
    condensed = []
    for gene, group in merged.groupby("gene_name"):
        top_row = group.iloc[0]
        top_cluster = (
            group.sort_values("mean_expression", ascending=False).iloc[0][
                ["cluster", "mean_expression"]
            ]
            if group["mean_expression"].notna().any()
            else {"cluster": None, "mean_expression": None}
        )
        condensed.append(
            {
                "gene": gene,
                "gene_id": top_row["gene_id"],
                "log2fc": float(top_row["log2fc"]),
                "padj": float(top_row["padj"]),
                "n_chip_peaks": int(top_row["n_peaks"]),
                "max_chip_score": float(top_row["max_score"]),
                "scrna_top_cluster": top_cluster.get("cluster"),
                "scrna_top_cluster_mean": (
                    float(top_cluster.get("mean_expression"))
                    if top_cluster.get("mean_expression") is not None
                    else None
                ),
            }
        )

    summary_json = {
        "n_genes_integrated": len(condensed),
        "top_hits": condensed[: args.top_n],
        "figures": {
            "heatmap": str(heatmap_path.name),
            "chip_bar": str(bar_path.name),
        },
    }
    Path(args.out_json).write_text(json.dumps(summary_json, indent=2))


def _plot_heatmap(merged: pd.DataFrame, path: Path) -> None:
    top_genes = merged.sort_values("padj").head(15)["gene_name"].unique()
    subset = merged[merged["gene_name"].isin(top_genes)]
    if subset.empty:
        plt.figure(figsize=(6, 3))
        plt.text(0.5, 0.5, "No overlapping genes", ha="center", va="center")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(path, dpi=150)
        plt.close()
        return
    matrix = (
        subset.pivot_table(index="gene_name", columns="cluster", values="mean_expression", aggfunc="max")
        .fillna(0.0)
    )
    plt.figure(figsize=(10, max(4, len(matrix) * 0.4)))
    sns.heatmap(matrix, cmap="viridis")
    plt.title("scRNA-seq expression of STAT1-bound DEGs")
    plt.xlabel("PBMC cluster")
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()


def _plot_bar(merged: pd.DataFrame, path: Path) -> None:
    top = (
        merged.sort_values("padj")
        .drop_duplicates("gene_name")
        .sort_values("n_peaks", ascending=False)
        .head(15)
    )
    if top.empty:
        plt.figure(figsize=(6, 3))
        plt.text(0.5, 0.5, "No peaks available", ha="center", va="center")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(path, dpi=150)
        plt.close()
        return
    plt.figure(figsize=(10, 4))
    sns.barplot(data=top, x="gene_name", y="n_peaks", color="#4477AA")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Number of STAT1 peaks")
    plt.xlabel("Gene")
    plt.title("ChIP-seq support for top DEGs")
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()


if __name__ == "__main__":
    main()

