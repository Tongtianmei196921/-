from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import requests


def _post_mygene_symbols(symbols: list[str]) -> dict[str, str]:
    if not symbols:
        return {}
    response = requests.post(
        "https://mygene.info/v3/query",
        data={
            "q": ",".join(symbols),
            "scopes": "symbol",
            "fields": "entrezgene,symbol",
            "species": "human",
            "size": 1,
        },
        timeout=120,
    )
    response.raise_for_status()
    mapping: dict[str, str] = {}
    for item in response.json():
        query = str(item.get("query", "")).upper()
        entrez = item.get("entrezgene")
        if query and entrez is not None:
            mapping[query] = str(entrez)
    return mapping


def map_symbols_to_entrez(symbols: list[str], batch_size: int = 900) -> dict[str, str]:
    mapping: dict[str, str] = {}
    unique = list(dict.fromkeys(str(symbol).upper() for symbol in symbols if str(symbol).strip()))
    for start in range(0, len(unique), batch_size):
        mapping.update(_post_mygene_symbols(unique[start : start + batch_size]))
    return mapping


def _write_gmt(path: Path, rows: dict[str, list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for name, genes in rows.items():
            handle.write("\t".join([name, "DrugReflector signature export", *genes]) + "\n")


def export_query_sets(signature_path: Path, out_dir: Path, top_n: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    signature = pd.read_csv(signature_path, index_col=0)

    up_symbols: dict[str, list[str]] = {}
    down_symbols: dict[str, list[str]] = {}
    all_symbols: list[str] = []
    for query_name, row in signature.iterrows():
        ordered = row.dropna().sort_values(ascending=False)
        up = ordered[ordered > 0].head(top_n).index.astype(str).str.upper().tolist()
        down = ordered[ordered < 0].sort_values(ascending=True).head(top_n).index.astype(str).str.upper().tolist()
        up_symbols[str(query_name)] = up
        down_symbols[str(query_name)] = down
        all_symbols.extend(up)
        all_symbols.extend(down)

    _write_gmt(out_dir / "clue_uptag_symbols.gmt", up_symbols)
    _write_gmt(out_dir / "clue_dntag_symbols.gmt", down_symbols)

    mapping = map_symbols_to_entrez(all_symbols)
    up_entrez = {
        query: [mapping[gene] for gene in genes if gene in mapping]
        for query, genes in up_symbols.items()
    }
    down_entrez = {
        query: [mapping[gene] for gene in genes if gene in mapping]
        for query, genes in down_symbols.items()
    }
    _write_gmt(out_dir / "clue_uptag_entrez.gmt", up_entrez)
    _write_gmt(out_dir / "clue_dntag_entrez.gmt", down_entrez)

    stats = []
    for query in up_symbols:
        stats.append(
            {
                "query": query,
                "up_symbols": len(up_symbols[query]),
                "down_symbols": len(down_symbols[query]),
                "up_entrez_mapped": len(up_entrez[query]),
                "down_entrez_mapped": len(down_entrez[query]),
            }
        )
    pd.DataFrame(stats).to_csv(out_dir / "clue_query_gene_set_mapping_stats.csv", index=False)
    (out_dir / "symbol_to_entrez_mapping.json").write_text(
        json.dumps(mapping, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export DrugReflector signature rows as CLUE/L1000 up/down GMT query files."
    )
    parser.add_argument(
        "--signature",
        type=Path,
        default=Path("outputs/GSE198138_processed/differential/GSE198138_DrugReflector_signature_matrix.csv"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("outputs/GSE198138_processed/clue_query"),
    )
    parser.add_argument("--top-n", type=int, default=150)
    args = parser.parse_args()
    export_query_sets(args.signature, args.out_dir, args.top_n)
    print(args.out_dir.resolve())


if __name__ == "__main__":
    main()
