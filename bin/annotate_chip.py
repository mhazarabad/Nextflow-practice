import argparse
import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd


def main(argv = None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--peaks", required=True, help="Path to BED (narrowPeak) file.")
    parser.add_argument("--gtf", required=True, help="Path to gene annotation GTF.gz.")
    parser.add_argument("--max-peaks", type=int, default=1000)
    parser.add_argument("--out-annotated", required=True)
    parser.add_argument("--out-genes", required=True)
    args = parser.parse_args(list(argv) if argv is not None else None)

    peaks_df = _load_peaks(Path(args.peaks), args.max_peaks)
    gene_intervals = _load_genes(Path(args.gtf))
    annotated = _annotate(peaks_df, gene_intervals)
    annotated.to_csv(args.out_annotated, sep="\t", index=False)

    agg = (
        annotated.groupby(["gene_id", "gene_name"])
        .agg(
            n_peaks=("peak_id", "count"),
            max_score=("score", "max"),
            mean_score=("score", "mean"),
        )
        .reset_index()
        .sort_values("n_peaks", ascending=False)
    )
    agg.to_csv(args.out_genes, sep="\t", index=False)

def _load_peaks(path: Path, max_peaks: int):
    cols = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "signal",
        "pvalue",
        "qvalue",
        "peak",
    ]
    if path.suffix == ".gz":
        opener = lambda p: gzip.open(p, "rt")
    else:
        opener = lambda p: open(p, "r")
    with opener(path) as handle:
        df = pd.read_csv(handle, sep="\t", names=cols, comment="#")
    df = df.sort_values("score", ascending=False).head(max_peaks).copy()
    df["peak_id"] = [f"peak_{i:04d}" for i in range(1, len(df) + 1)]
    return df[["peak_id", "chrom", "start", "end", "score"]]

def _load_genes(gtf_path: Path):
    intervals: Dict[str, List[Tuple[int, int, str, str]]] = defaultdict(list)
    opener = gzip.open if gtf_path.suffix == ".gz" else open
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            attrs = _parse_gtf_attributes(parts[8])
            gene_id = attrs.get("gene_id")
            gene_name = attrs.get("gene_name", gene_id)
            if gene_id is None:
                continue
            intervals[chrom].append((start, end, gene_id.split(".")[0], gene_name))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals

def _parse_gtf_attributes(field: str):
    attrs: Dict[str, str] = {}
    for item in field.strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        attrs[key] = value.strip('"')
    return attrs

def _annotate(
    peaks: pd.DataFrame, gene_intervals: Dict[str, List[Tuple[int, int, str, str]]]
):
    records: List[Dict[str, object]] = []
    for _, peak in peaks.iterrows():
        chrom = peak["chrom"]
        if chrom not in gene_intervals:
            continue
        start, end = int(peak["start"]), int(peak["end"])
        for gene_start, gene_end, gene_id, gene_name in gene_intervals[chrom]:
            if gene_end < start:
                continue
            if gene_start > end:
                break
            if gene_end >= start and gene_start <= end:
                records.append(
                    {
                        "peak_id": peak["peak_id"],
                        "chrom": chrom,
                        "peak_start": start,
                        "peak_end": end,
                        "score": peak["score"],
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                    }
                )
    return pd.DataFrame(records)

main()

