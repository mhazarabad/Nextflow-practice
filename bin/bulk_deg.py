import argparse
from pathlib import Path
from typing import Dict
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection


def main(argv = None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", required=True, help="TSV gene x sample matrix.")
    parser.add_argument("--metadata", required=True, help="Sample metadata TSV.")
    parser.add_argument("--gtf", required=True, help="Gene annotation GTF.gz for mapping IDs.")
    parser.add_argument("--out-deg", required=True)
    parser.add_argument("--out-toplist", required=True)
    parser.add_argument("--min-logfc", type=float, default=0.5)
    parser.add_argument("--padj-threshold", type=float, default=0.1)
    args = parser.parse_args(list(argv) if argv is not None else None)

    expr = pd.read_csv(args.expression, sep="\t", index_col=0)
    meta = pd.read_csv(args.metadata, sep="\t")
    if "condition" not in meta.columns:
        raise RuntimeError("Metadata missing 'condition' column.")

    groups = meta["condition"].unique()
    if len(groups) != 2:
        raise RuntimeError(f"Expected exactly two conditions, found {groups}.")

    group_a, group_b = groups
    samples_a = meta.loc[meta["condition"] == group_a, "sample_id"].tolist()
    samples_b = meta.loc[meta["condition"] == group_b, "sample_id"].tolist()
    if len(samples_a) < 2 or len(samples_b) < 2:
        raise RuntimeError("Need at least two replicates per condition.")

    expr_a = expr[samples_a].astype(float)
    expr_b = expr[samples_b].astype(float)

    mean_a = expr_a.mean(axis=1)
    mean_b = expr_b.mean(axis=1)
    log2fc = np.log2((mean_b + 1e-3) / (mean_a + 1e-3))

    t_stat, pvals = stats.ttest_ind(expr_b, expr_a, axis=1, equal_var=False)
    pvals = np.nan_to_num(pvals, nan=1.0)
    _, padj = fdrcorrection(pvals, alpha=args.padj_threshold)

    deg = pd.DataFrame(
        {
            "gene_id": expr.index,
            "mean_a": mean_a.values,
            "mean_b": mean_b.values,
            "log2fc": log2fc.values,
            "pvalue": pvals,
            "padj": padj,
        }
    )
    deg = _annotate_gene_names(deg, Path(args.gtf))
    deg = deg.sort_values("padj").reset_index(drop=True)
    deg.to_csv(args.out_deg, sep="\t", index=False)

    toplist = deg[
        (deg["padj"] <= args.padj_threshold) & (deg["log2fc"].abs() >= args.min_logfc)
    ].copy()
    toplist.to_csv(args.out_toplist, sep="\t", index=False)

def _annotate_gene_names(deg: pd.DataFrame, gtf_path: Path):
    mapping = {}
    with gzip_open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = _parse_gtf_attributes(parts[8])
            gene_id = attrs.get("gene_id")
            gene_name = attrs.get("gene_name")
            if gene_id and gene_name:
                mapping[gene_id.split(".")[0]] = gene_name
    deg["gene_id_stripped"] = deg["gene_id"].str.split(".").str[0]
    deg["gene_name"] = deg["gene_id_stripped"].map(mapping)
    deg.drop(columns=["gene_id_stripped"], inplace=True)
    return deg

def _parse_gtf_attributes(field: str):
    entries = {}
    for item in field.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if " " not in item:
            continue
        key, value = item.split(" ", 1)
        entries[key] = value.strip('"')
    return entries

def gzip_open(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        import gzip

        return gzip.open(path, mode)
    return path.open(mode)


main()