from __future__ import annotations

import gzip
import io
import re
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import requests
from scipy import stats

from drugreflector.drug_reflector import DrugReflector
from drugreflector.repo_annotations import get_compound_annotations
from drugreflector.utils import compute_vscores_adata


ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / ".geo_cache"
OUT = ROOT / "outputs" / "GSE198138_processed"
RAW_TAR = CACHE / "GSE198138_RAW.tar"
SERIES_MATRIX = CACHE / "GSE198138_series_matrix.txt.gz"
ENSEMBL_MAP = CACHE / "ensembl_human_gene_symbol.tsv"


def _clean_geo_value(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        value = value[1:-1]
    return value


def _read_series_metadata() -> pd.DataFrame:
    records: list[dict[str, str]] = []
    repeated_characteristics: list[list[str]] = []

    with gzip.open(SERIES_MATRIX, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            parts = [_clean_geo_value(part) for part in line.rstrip("\n").split("\t")]
            key = parts[0]
            values = parts[1:]
            if key == "!Sample_geo_accession":
                records = [{"geo_accession": value} for value in values]
            elif key == "!Sample_title" and records:
                for record, value in zip(records, values):
                    record["title"] = value
            elif key == "!Sample_source_name_ch1" and records:
                for record, value in zip(records, values):
                    record["source_name"] = value
            elif key == "!Sample_characteristics_ch1" and records:
                repeated_characteristics.append(values)

    for values in repeated_characteristics:
        for record, value in zip(records, values):
            if ":" not in value:
                continue
            field, raw = value.split(":", 1)
            record[field.strip().lower().replace(" ", "_")] = raw.strip()

    metadata = pd.DataFrame(records)
    if metadata.empty:
        raise ValueError("No sample metadata could be recovered from the series matrix.")
    metadata["disease_status"] = metadata["disease_status"].replace(
        {"Unaffected": "Unaffected", "FXS": "FXS"}
    )
    metadata["differentiation_day"] = metadata["differentiation_day"].str.replace(" ", "_", regex=False)
    return metadata.set_index("geo_accession", drop=False)


def _read_inner_csv(inner: tarfile.TarFile, suffix: str) -> pd.DataFrame:
    matches = [member for member in inner.getmembers() if member.name.endswith(suffix)]
    if not matches:
        raise ValueError(f"Missing {suffix} in {inner.name}")
    with inner.extractfile(matches[0]) as handle:
        if handle is None:
            raise ValueError(f"Could not open {matches[0].name}")
        return pd.read_csv(handle)


def _sample_expression_from_analysis(raw_bytes: bytes) -> tuple[pd.Series, int, int]:
    with tarfile.open(fileobj=io.BytesIO(raw_bytes), mode="r:gz") as inner:
        clusters = _read_inner_csv(inner, "clustering/graphclust/clusters.csv")
        diffexp = _read_inner_csv(inner, "diffexp/graphclust/differential_expression.csv")

    if "Cluster" not in clusters.columns:
        raise ValueError("Cell Ranger graphclust clusters.csv has no Cluster column.")
    if "Gene ID" not in diffexp.columns:
        raise ValueError("Cell Ranger differential_expression.csv has no Gene ID column.")

    cluster_counts = clusters["Cluster"].astype(str).value_counts()
    total_cells = int(cluster_counts.sum())
    mean_columns = []
    for column in diffexp.columns:
        match = re.match(r"Cluster\s+(\d+)\s+Mean UMI Counts$", str(column))
        if match:
            mean_columns.append((match.group(1), column))
    if not mean_columns:
        raise ValueError("No cluster mean UMI columns were found.")

    weighted = np.zeros(diffexp.shape[0], dtype=float)
    used_clusters = 0
    for cluster_id, column in mean_columns:
        weight = int(cluster_counts.get(cluster_id, 0))
        if weight <= 0:
            continue
        weighted += pd.to_numeric(diffexp[column], errors="coerce").fillna(0).to_numpy(dtype=float) * weight
        used_clusters += 1

    weighted = weighted / max(total_cells, 1)
    series = pd.Series(weighted, index=diffexp["Gene ID"].astype(str).str.replace(r"\.\d+$", "", regex=True))
    series = series.groupby(level=0).mean()
    return series, total_cells, used_clusters


def build_sample_expression() -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = _read_series_metadata()
    series_by_sample: dict[str, pd.Series] = {}
    cell_counts: dict[str, int] = {}
    used_clusters: dict[str, int] = {}

    with tarfile.open(RAW_TAR, mode="r") as outer:
        for member in outer.getmembers():
            if not member.name.endswith("_analysis.tar.gz"):
                continue
            accession_match = re.match(r"(GSM\d+)_", member.name)
            if not accession_match:
                continue
            accession = accession_match.group(1)
            with outer.extractfile(member) as handle:
                if handle is None:
                    raise ValueError(f"Could not open {member.name}")
                expression, n_cells, n_clusters = _sample_expression_from_analysis(handle.read())
            series_by_sample[accession] = expression
            cell_counts[accession] = n_cells
            used_clusters[accession] = n_clusters

    expression = pd.DataFrame(series_by_sample).sort_index()
    sample_meta = metadata.loc[expression.columns].copy()
    sample_meta["n_cells_graphclust"] = pd.Series(cell_counts)
    sample_meta["n_graphclust_clusters_used"] = pd.Series(used_clusters)
    return expression, sample_meta


def download_ensembl_mapping() -> pd.DataFrame:
    if not ENSEMBL_MAP.exists():
        query = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id" />
    <Attribute name="hgnc_symbol" />
    <Attribute name="external_gene_name" />
  </Dataset>
</Query>"""
        response = requests.post(
            "https://www.ensembl.org/biomart/martservice",
            data={"query": query},
            timeout=120,
        )
        response.raise_for_status()
        ENSEMBL_MAP.write_text(response.text, encoding="utf-8")

    mapping = pd.read_csv(ENSEMBL_MAP, sep="\t")
    mapping.columns = ["ensembl_gene_id", "hgnc_symbol", "external_gene_name"]
    mapping["ensembl_gene_id"] = mapping["ensembl_gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    mapping["symbol"] = mapping["hgnc_symbol"].fillna("").astype(str).str.strip()
    fallback = mapping["external_gene_name"].fillna("").astype(str).str.strip()
    mapping.loc[mapping["symbol"].eq("") & ~fallback.str.startswith("ENSG"), "symbol"] = fallback
    mapping = mapping[mapping["symbol"].ne("")].drop_duplicates("ensembl_gene_id")
    return mapping[["ensembl_gene_id", "symbol"]]


def map_to_hgnc(expression: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int | float]]:
    mapping = download_ensembl_mapping()
    mapped = expression.merge(mapping, left_index=True, right_on="ensembl_gene_id", how="inner")
    collapsed = mapped.drop(columns=["ensembl_gene_id"]).groupby("symbol").mean()
    stats_dict = {
        "input_ensembl_genes": int(expression.shape[0]),
        "mapped_gene_rows": int(mapped["ensembl_gene_id"].nunique()),
        "unique_hgnc_symbols": int(collapsed.shape[0]),
        "mapping_rate_percent": round(float(mapped["ensembl_gene_id"].nunique() / max(expression.shape[0], 1) * 100), 2),
    }
    return collapsed.sort_index(), stats_dict


def _bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    out = np.full_like(pvalues, np.nan, dtype=float)
    valid = np.isfinite(pvalues)
    if not valid.any():
        return out
    p = pvalues[valid]
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / (np.arange(len(ranked)) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    valid_idx = np.flatnonzero(valid)
    out[valid_idx[order]] = adjusted
    return out


def differential_table(expression_symbols: pd.DataFrame, sample_meta: pd.DataFrame, label: str, mask: pd.Series) -> pd.DataFrame:
    samples = sample_meta.index[mask].tolist()
    subset_expr = expression_symbols[samples].T
    subset_meta = sample_meta.loc[samples]

    adata = ad.AnnData(
        X=subset_expr.to_numpy(dtype=float),
        obs=subset_meta.copy(),
        var=pd.DataFrame(index=subset_expr.columns),
    )
    vscores = compute_vscores_adata(adata, "disease_status", "Unaffected", "FXS")

    control = subset_expr.loc[subset_meta["disease_status"] == "Unaffected"]
    case = subset_expr.loc[subset_meta["disease_status"] == "FXS"]
    pseudo = max(float(np.nanmedian(subset_expr.to_numpy()) * 0.01), 1e-6)

    mean_control = control.mean(axis=0)
    mean_case = case.mean(axis=0)
    log2fc = np.log2((mean_case + pseudo) / (mean_control + pseudo))
    ttest = stats.ttest_ind(case, control, axis=0, equal_var=False, nan_policy="omit")
    pvalue = np.asarray(ttest.pvalue, dtype=float)

    table = pd.DataFrame(
        {
            "gene": subset_expr.columns,
            "v_score_FXS_vs_Unaffected": vscores.reindex(subset_expr.columns).to_numpy(),
            "log2FC_FXS_vs_Unaffected": log2fc.reindex(subset_expr.columns).to_numpy(),
            "mean_FXS": mean_case.reindex(subset_expr.columns).to_numpy(),
            "mean_Unaffected": mean_control.reindex(subset_expr.columns).to_numpy(),
            "pvalue_welch_on_derived_sample_means": pvalue,
            "padj_BH_exploratory": _bh_adjust(pvalue),
            "n_FXS_samples": int(case.shape[0]),
            "n_Unaffected_samples": int(control.shape[0]),
            "contrast": label,
        }
    )
    table = table.sort_values("v_score_FXS_vs_Unaffected", ascending=False)
    return table


def flatten_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    annotations = get_compound_annotations(predictions.index.astype(str).tolist())
    for compound in predictions.index.astype(str):
        annotation = annotations.get(compound, {})
        for signature in predictions.columns.get_level_values(1).unique():
            rows.append(
                {
                    "signature": signature,
                    "compound": compound,
                    "display_name": annotation.get("display_name"),
                    "chemical_name": annotation.get("chemical_name"),
                    "target": annotation.get("target"),
                    "moa": annotation.get("moa"),
                    "rank": int(predictions.loc[compound, ("rank", signature)]) + 1,
                    "logit": float(predictions.loc[compound, ("logit", signature)]),
                    "probability": float(predictions.loc[compound, ("prob", signature)]),
                    "pubchem_cid": annotation.get("pubchem_cid"),
                    "pubchem_url": annotation.get("pubchem_url"),
                }
            )
    return pd.DataFrame(rows).sort_values(["signature", "rank"])


def write_markdown_report(
    sample_meta: pd.DataFrame,
    mapping_stats: dict[str, int | float],
    differential_tables: dict[str, pd.DataFrame],
    predictions_long: pd.DataFrame,
) -> None:
    lines = [
        "# GSE198138 DrugReflector 探索性处理报告",
        "",
        "## 数据边界",
        "",
        "- GEO：GSE198138",
        "- 物种：Homo sapiens",
        "- 原始公开文件：GSE198138_RAW.tar",
        "- 可用内容：Cell Ranger `analysis.tar.gz` 与 `.cloupe.gz`",
        "- 未发现标准 10x count matrix：未发现 `matrix.mtx`、`features.tsv`、`barcodes.tsv`、`filtered_feature_bc_matrix.h5` 或 `.h5ad`",
        "- 本次处理方式：由每个样本的 graphclust cluster mean UMI 与 cluster cell counts 重建样本级平均表达，再计算 FXS 相对 Unaffected 的探索性 signature",
        "- 严格限制：这不是原始单细胞 counts 的标准重分析，也不能替代作者原始 Seurat/Scanpy 流程或实验验证",
        "",
        "## 样本概况",
        "",
        sample_meta[
            [
                column
                for column in [
                    "title",
                    "geo_accession",
                    "source_name",
                    "disease_status",
                    "differentiation_day",
                    "n_cells_graphclust",
                ]
                if column in sample_meta.columns
            ]
        ].to_markdown(),
        "",
        "## 基因映射",
        "",
        f"- 输入 Ensembl genes：{mapping_stats['input_ensembl_genes']}",
        f"- 成功映射 gene rows：{mapping_stats['mapped_gene_rows']}",
        f"- 唯一 HGNC symbols：{mapping_stats['unique_hgnc_symbols']}",
        f"- 映射率：{mapping_stats['mapping_rate_percent']}%",
        "",
    ]

    for label, table in differential_tables.items():
        safe = label.replace(" ", "_")
        top_up = table.sort_values("v_score_FXS_vs_Unaffected", ascending=False).head(15)
        top_down = table.sort_values("v_score_FXS_vs_Unaffected", ascending=True).head(15)
        lines.extend(
            [
                f"## {label} 差异 signature",
                "",
                "方向定义：正值表示 FXS 高于 Unaffected；负值表示 FXS 低于 Unaffected。",
                "",
                "### Top 15 FXS-up candidates",
                "",
                top_up[["gene", "v_score_FXS_vs_Unaffected", "log2FC_FXS_vs_Unaffected", "padj_BH_exploratory"]].to_markdown(index=False),
                "",
                "### Top 15 FXS-down candidates",
                "",
                top_down[["gene", "v_score_FXS_vs_Unaffected", "log2FC_FXS_vs_Unaffected", "padj_BH_exploratory"]].to_markdown(index=False),
                "",
                "### Top 15 DrugReflector compounds",
                "",
                predictions_long[predictions_long["signature"].eq(safe)]
                .head(15)[["rank", "compound", "display_name", "chemical_name", "probability", "target", "moa", "pubchem_url"]]
                .to_markdown(index=False),
                "",
            ]
        )

    lines.extend(
        [
            "## 下一步建议",
            "",
            "- 优先联系作者或从 SRA FASTQ 重新跑 Cell Ranger/STARsolo，以获得真正的 cell-by-gene count matrix。",
            "- 若继续使用本次探索性结果，应把它定位为候选药物初筛，不应作为正式论文中的最终差异分析证据。",
            "- 正式机制验证建议优先围绕 FXS vs Unaffected 的复现差异基因、DrugReflector 高排位化合物及其已知 target/MOA 的交集设计实验。",
        ]
    )
    (OUT / "GSE198138_processing_report.zh-CN.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    for subdir in ["expression", "differential", "predictions"]:
        (OUT / subdir).mkdir(parents=True, exist_ok=True)

    expression_ensembl, sample_meta = build_sample_expression()
    expression_symbols, mapping_stats = map_to_hgnc(expression_ensembl)

    expression_symbols.T.to_csv(OUT / "expression" / "GSE198138_derived_sample_expression_HGNC.csv")
    sample_meta.to_csv(OUT / "expression" / "GSE198138_sample_metadata.csv")

    expression_adata = ad.AnnData(
        X=expression_symbols.T.to_numpy(dtype=float),
        obs=sample_meta.copy(),
        var=pd.DataFrame(index=expression_symbols.index),
    )
    expression_adata.uns["source_note"] = (
        "Derived from Cell Ranger graphclust cluster mean UMI counts weighted by graphclust cluster cell counts; "
        "not a raw single-cell count matrix."
    )
    expression_adata.write_h5ad(OUT / "expression" / "GSE198138_derived_sample_expression_HGNC.h5ad", compression="gzip")

    contrasts = {
        "FXS_vs_Unaffected_all": pd.Series(True, index=sample_meta.index),
        "FXS_vs_Unaffected_day22": sample_meta["differentiation_day"].eq("day_22"),
        "FXS_vs_Unaffected_day42_48": sample_meta["differentiation_day"].eq("day_42-48"),
    }
    differential_tables: dict[str, pd.DataFrame] = {}
    signature_rows: dict[str, pd.Series] = {}
    for label, mask in contrasts.items():
        table = differential_table(expression_symbols, sample_meta, label, mask)
        differential_tables[label] = table
        table.to_csv(OUT / "differential" / f"GSE198138_{label}_derived_differential_table.csv", index=False)
        signature_rows[label] = table.set_index("gene")["v_score_FXS_vs_Unaffected"].sort_index()

    signature_matrix = pd.DataFrame(signature_rows).T.fillna(0)
    signature_matrix.to_csv(OUT / "differential" / "GSE198138_DrugReflector_signature_matrix.csv")
    signature_adata = ad.AnnData(
        X=signature_matrix.to_numpy(dtype=float),
        obs=pd.DataFrame(index=signature_matrix.index),
        var=pd.DataFrame(index=signature_matrix.columns),
    )
    signature_adata.write_h5ad(OUT / "differential" / "GSE198138_DrugReflector_signature_matrix.h5ad", compression="gzip")

    checkpoints = [ROOT / "checkpoints" / f"model_fold_{idx}.pt" for idx in range(3)]
    model = DrugReflector([str(path) for path in checkpoints])
    predictions = model.predict(signature_matrix, n_top=50)
    predictions.to_csv(OUT / "predictions" / "GSE198138_DrugReflector_predictions_raw.csv")
    predictions_long = flatten_predictions(predictions)
    predictions_long.to_csv(OUT / "predictions" / "GSE198138_DrugReflector_predictions_top50_annotated.csv", index=False)

    write_markdown_report(sample_meta, mapping_stats, differential_tables, predictions_long)
    print(f"Wrote processed GSE198138 outputs to {OUT}")


if __name__ == "__main__":
    main()
