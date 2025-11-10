import argparse
from pathlib import Path
from typing import Iterable
import scanpy as sc


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-file", required=True)
    args = parser.parse_args(list(argv) if argv is not None else None)

    adata = sc.datasets.pbmc3k_processed()
    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)


main()

