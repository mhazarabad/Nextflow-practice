#!/usr/bin/env python3
"""
Download GEO RNA-seq normalized counts and sample metadata, then subset to a
specific immune cell type and disease contrast.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd
import requests


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        return dest
    with requests.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        with dest.open("wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 16):
                if chunk:
                    fh.write(chunk)
    return dest


def _parse_series_matrix(path: Path) -> Dict[str, Dict[str, str]]:
    sample_titles: List[str] = []
    sample_gsms: List[str] = []
    characteristic_rows: List[List[str]] = []

    line_pattern = re.compile(r'^"!?(.*?)"$')

    with gzip.open(path, "rt") as handle:
        for raw in handle:
            if raw.startswith("!Sample_title"):
                sample_titles = _split_values(raw)
            elif raw.startswith("!Sample_geo_accession"):
                sample_gsms = _split_values(raw)
            elif raw.startswith("!Sample_characteristics_ch1"):
                characteristic_rows.append(_split_values(raw))

    if not sample_titles:
        raise RuntimeError("Failed to parse sample titles from series matrix.")
    if len(sample_titles) != len(sample_gsms):
        raise RuntimeError("Mismatch between titles and GSM identifiers.")

    samples = [defaultdict(str) for _ in sample_titles]
    for idx, gsm in enumerate(sample_gsms):
        samples[idx]["gsm"] = gsm
        samples[idx]["sample_id"] = sample_titles[idx]

    for row in characteristic_rows:
        for idx, value in enumerate(row):
            if value in {"--", "", "NA"}:
                continue
            if ": " in value:
                key, val = value.split(": ", 1)
            elif ":" in value:
                key, val = value.split(":", 1)
            else:
                key, val = f"attr_{len(samples[idx])}", value
            key = key.strip().lower().replace(" ", "_")
            samples[idx][key] = val.strip()

    metadata: Dict[str, Dict[str, str]] = {}
    for sample in samples:
        metadata[sample["sample_id"]] = sample
    return metadata


def _split_values(line: str) -> List[str]:
    parts = line.strip().split("\t")[1:]
    cleaned: List[str] = []
    for entry in parts:
        entry = entry.strip().strip('"')
        cleaned.append(entry)
    return cleaned


def _load_expression(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        df = pd.read_csv(handle, sep="\t")
    if "genenames" not in df.columns:
        raise RuntimeError("Counts table is missing 'genenames' column.")
    df = df.rename(columns={"genenames": "gene_id"})
    df = df.set_index("gene_id")
    return df


def _filter_expression(
    expression: pd.DataFrame,
    metadata: Dict[str, Dict[str, str]],
    celltype: str,
    include_conditions: Tuple[str, str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    harmonized = {}
    for sample_id, attrs in metadata.items():
        if attrs.get("celltype", "").lower() != celltype.lower():
            continue
        condition = attrs.get("diseasestatus") or attrs.get("disease_status")
        if condition is None:
            continue
        if condition not in include_conditions:
            continue
        attrs = dict(attrs)
        attrs["condition"] = condition
        harmonized[sample_id] = attrs

    if len(harmonized) == 0:
        raise RuntimeError(
            f"No samples matched celltype='{celltype}' with conditions={include_conditions}."
        )

    missing = [sid for sid in harmonized if sid not in expression.columns]
    if missing:
        raise RuntimeError(f"Expression table missing samples: {missing[:5]}")

    filtered_expression = expression[list(harmonized.keys())].copy()
    metadata_df = pd.DataFrame.from_dict(harmonized, orient="index")
    # drop non-informative columns (all identical)
    metadata_df = metadata_df.loc[:, metadata_df.nunique(dropna=False) > 1].copy()
    metadata_df["sample_id"] = metadata_df.index
    cols = ["sample_id", "gsm", "celltype", "condition", "donorid"]
    metadata_df = metadata_df.reindex(
        columns=cols + [c for c in metadata_df.columns if c not in cols]
    )
    return filtered_expression, metadata_df


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix-url", required=True)
    parser.add_argument("--counts-url", required=True)
    parser.add_argument("--celltype", required=True)
    parser.add_argument(
        "--conditions",
        required=True,
        help="Comma-separated list of two condition labels to keep.",
    )
    parser.add_argument("--out-expression", required=True)
    parser.add_argument("--out-metadata", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--cache-dir", default="cache")
    args = parser.parse_args(list(argv) if argv is not None else None)

    cache_dir = Path(args.cache_dir)
    matrix_path = _download(args.matrix_url, cache_dir / "series_matrix.txt.gz")
    counts_path = _download(args.counts_url, cache_dir / "normalized_counts.txt.gz")

    metadata = _parse_series_matrix(matrix_path)
    expression = _load_expression(counts_path)

    condition_pairs = tuple(item.strip() for item in args.conditions.split(",") if item.strip())
    if len(condition_pairs) != 2:
        raise ValueError("Exactly two conditions must be provided.")

    filtered_expression, metadata_df = _filter_expression(
        expression, metadata, args.celltype, condition_pairs  # type: ignore[arg-type]
    )

    filtered_expression = filtered_expression.loc[
        (filtered_expression > 0).sum(axis=1) > 0
    ]

    filtered_expression.to_csv(args.out_expression, sep="\t")
    metadata_df.to_csv(args.out_metadata, sep="\t", index=False)

    summary = {
        "celltype": args.celltype,
        "conditions": list(condition_pairs),
        "n_samples": int(metadata_df.shape[0]),
        "n_genes": int(filtered_expression.shape[0]),
        "samples_per_condition": metadata_df.groupby("condition")["sample_id"].count().to_dict(),
    }
    Path(args.summary_json).write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

