import argparse
from pathlib import Path
from typing import Iterable
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def main(argv = None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--deg", required=True)
    parser.add_argument("--chip-genes", required=True)
    parser.add_argument("--scrna-signature", required=True)
    parser.add_argument("--out-summary", required=True)
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

def _plot_heatmap(merged: pd.DataFrame, path: Path):
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

def _plot_bar(merged: pd.DataFrame, path: Path):
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




main()
