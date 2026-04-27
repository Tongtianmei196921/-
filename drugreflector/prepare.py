"""Objective input inspection and preparation for DrugReflector."""

from __future__ import annotations

import io
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData

from .utils import compute_vscores_adata, load_h5ad_file, pseudobulk_adata

GENE_COLUMN_HINTS = (
    "gene",
    "gene_symbol",
    "genesymbol",
    "symbol",
    "hgnc",
    "gene_name",
    "external_gene_name",
)
SCORE_COLUMN_HINTS = (
    "vscore",
    "v_score",
    "score",
    "log2foldchange",
    "logfoldchange",
    "logfc",
    "avg_log2fc",
    "avg_logfc",
    "stat",
    "t",
    "waldstat",
    "coef",
)
EXCLUDE_SCORE_HINTS = ("pvalue", "padj", "qvalue", "fdr", "pct", "percent")
GROUP_COLUMN_HINTS = (
    "group",
    "condition",
    "cell_type",
    "celltype",
    "cluster",
    "class",
    "status",
    "phenotype",
    "disease",
    "treatment",
)
CONTROL_HINTS = (
    "control",
    "normal",
    "healthy",
    "untreated",
    "vehicle",
    "wildtype",
    "wild type",
    "parental",
    "sensitive",
    "baseline",
    "adjacent",
    "benign",
    "non-tumor",
    "nontumor",
    "mock",
)
CASE_HINTS = (
    "tumor",
    "cancer",
    "disease",
    "treated",
    "resistant",
    "mutant",
    "knockout",
    "ko",
    "kd",
    "patient",
    "infected",
    "stimulated",
    "drug",
)
SAMPLE_ID_COLUMN_HINTS = (
    "sample",
    "sample_id",
    "sampleid",
    "donor",
    "patient",
    "individual",
    "replicate",
    "orig_ident",
    "orig.ident",
    "library",
    "library_id",
    "subject",
    "biosample",
)


@dataclass
class PreparedInput:
    adata: AnnData | None
    summary: dict[str, object]


def _canonical(text: str) -> str:
    return "".join(ch for ch in str(text).strip().lower() if ch.isalnum() or ch == "_")


def _normalize_gene_name(name: object) -> str:
    gene = str(name).strip().upper()
    gene = gene.split("///")[0].split("//")[0].split(";")[0].split(",")[0].strip()
    gene = gene.removesuffix("_AT")
    if "." in gene and gene.startswith("ENSG"):
        gene = gene.split(".", 1)[0]
    elif "." in gene:
        gene = gene.split(".", 1)[0]
    gene = "".join(ch for ch in gene if ch.isalnum() or ch == "-")
    return gene


def _looks_like_gene_name(value: object) -> bool:
    gene = _normalize_gene_name(value)
    if len(gene) < 2 or len(gene) > 20:
        return False
    if gene.isdigit():
        return False
    if gene.startswith(("SAMPLE", "CELL", "PATIENT", "CTRL", "CASE")):
        return False
    return all(ch.isalnum() or ch == "-" for ch in gene)


def _gene_likeness(labels: pd.Index) -> float:
    if len(labels) == 0:
        return 0.0
    subset = labels.astype(str)[: min(len(labels), 300)]
    return sum(_looks_like_gene_name(value) for value in subset) / len(subset)


def _is_raw_expression(values: np.ndarray) -> bool:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return False
    if np.nanmin(finite) < 0:
        return False
    if np.nanmax(finite) > 50:
        return True
    return np.nanpercentile(finite, 95) > 10


def _sample_label_grouping(labels: pd.Index) -> tuple[str, str] | None:
    if len(labels) < 4:
        return None

    tokens = []
    for raw in labels.astype(str):
        lower = raw.lower()
        token = (
            lower.replace("-", "_")
            .replace(" ", "_")
            .split("_")[0]
            .rstrip("0123456789")
        )
        if not token or token in {"sample", "cell", "x"}:
            return None
        tokens.append(token)

    counts = pd.Series(tokens).value_counts()
    if len(counts) != 2 or counts.min() < 2:
        return None

    group1, group2 = counts.index.tolist()
    joined = " ".join((group1, group2))
    if not any(hint in joined for hint in CONTROL_HINTS) and not any(
        hint in joined for hint in CASE_HINTS
    ):
        return None
    score1 = sum(hint in group1 for hint in CONTROL_HINTS) - sum(hint in group1 for hint in CASE_HINTS)
    score2 = sum(hint in group2 for hint in CONTROL_HINTS) - sum(hint in group2 for hint in CASE_HINTS)
    if score1 < score2:
        group1, group2 = group2, group1
    return group1, group2


def _build_group_series_from_labels(labels: pd.Index) -> pd.Series | None:
    grouping = _sample_label_grouping(labels)
    if not grouping:
        return None
    group1, group2 = grouping
    values = []
    for raw in labels.astype(str):
        lower = raw.lower().replace("-", "_").replace(" ", "_")
        token = lower.split("_")[0].rstrip("0123456789")
        values.append(token)
    return pd.Series(values, index=labels, name="_label_group")


def _choose_gene_column(columns: list[str]) -> str | None:
    normalized = {_canonical(col): col for col in columns}
    for hint in GENE_COLUMN_HINTS:
        if hint in normalized:
            return normalized[hint]
    for column in columns:
        canon = _canonical(column)
        if "gene" in canon or "symbol" in canon:
            return column
    return None


def _choose_score_column(frame: pd.DataFrame) -> str | None:
    numeric_columns = [col for col in frame.columns if pd.api.types.is_numeric_dtype(frame[col])]
    normalized = {_canonical(col): col for col in numeric_columns}
    for hint in SCORE_COLUMN_HINTS:
        if hint in normalized:
            return normalized[hint]
    filtered = [
        col
        for col in numeric_columns
        if not any(hint in _canonical(col) for hint in EXCLUDE_SCORE_HINTS)
    ]
    if len(filtered) == 1:
        return filtered[0]
    for col in filtered:
        canon = _canonical(col)
        if "log" in canon or "score" in canon or "stat" in canon:
            return col
    return None


def _series_to_adata(series: pd.Series, name: str) -> AnnData:
    clean = series.dropna().copy()
    clean.index = pd.Index([_normalize_gene_name(gene) for gene in clean.index])
    clean = clean[clean.index != ""]
    clean = clean.groupby(level=0).mean()
    return AnnData(
        X=clean.to_numpy(dtype=np.float32).reshape(1, -1),
        obs=pd.DataFrame(index=[name]),
        var=pd.DataFrame(index=clean.index),
    )


def _frame_to_adata(frame: pd.DataFrame) -> AnnData:
    clean = frame.copy()
    clean.columns = pd.Index([_normalize_gene_name(col) for col in clean.columns])
    clean = clean.loc[:, clean.columns != ""]
    clean = clean.T.groupby(level=0).mean().T
    return AnnData(
        X=clean.to_numpy(dtype=np.float32, copy=True),
        obs=pd.DataFrame(index=clean.index.astype(str)),
        var=pd.DataFrame(index=clean.columns.astype(str)),
    )


def _score_group_column(column: str, values: pd.Series) -> float:
    text = values.astype("string").fillna("").astype(str).str.strip()
    unique = [value for value in text.unique().tolist() if value]
    if len(unique) != 2:
        return -1.0

    score = 0.0
    canon = _canonical(column)
    has_column_hint = any(hint in canon for hint in GROUP_COLUMN_HINTS)
    if has_column_hint:
        score += 8.0
    joined = " ".join(value.lower() for value in unique)
    has_control_hint = any(hint in joined for hint in CONTROL_HINTS)
    has_case_hint = any(hint in joined for hint in CASE_HINTS)
    if has_control_hint:
        score += 6.0
    if has_case_hint:
        score += 6.0
    if not has_control_hint and not has_case_hint:
        return -1.0
    counts = text.value_counts()
    if counts.min() >= 2:
        score += 4.0
    return score


def _auto_group_from_obs(
    obs: pd.DataFrame,
    *,
    group_column: str | None = None,
    group1_value: str | None = None,
    group2_value: str | None = None,
) -> tuple[str, str, str] | None:
    if group_column and group1_value and group2_value:
        if group_column not in obs.columns:
            raise ValueError(f"Column '{group_column}' not found in metadata.")
        if group1_value == group2_value:
            raise ValueError("Selected comparison groups must be different.")
        text = obs[group_column].astype("string").fillna("").astype(str).str.strip()
        if not (text == group1_value).any():
            raise ValueError(
                f"No samples found with {group_column}='{group1_value}'."
            )
        if not (text == group2_value).any():
            raise ValueError(
                f"No samples found with {group_column}='{group2_value}'."
            )
        return group_column, group1_value, group2_value

    ranked: list[tuple[float, str]] = []
    for column in obs.columns:
        ranked.append((_score_group_column(column, obs[column]), column))
    ranked.sort(reverse=True)
    if not ranked or ranked[0][0] < 0:
        return None
    column = ranked[0][1]
    values = [value for value in obs[column].astype("string").fillna("").astype(str).unique().tolist() if value]
    if len(values) != 2:
        return None

    def _value_score(value: str) -> tuple[int, int, str]:
        lower = value.lower()
        return (
            0 if sum(hint in lower for hint in CONTROL_HINTS) >= sum(hint in lower for hint in CASE_HINTS) else 1,
            -sum(hint in lower for hint in CONTROL_HINTS),
            lower,
        )

    values = sorted(values, key=_value_score)
    return column, values[0], values[1]


def _candidate_group_columns(obs: pd.DataFrame) -> list[dict[str, object]]:
    candidates = []
    for column in obs.columns:
        text = obs[column].astype("string").fillna("").astype(str).str.strip()
        unique = [value for value in text.unique().tolist() if value]
        if 2 <= len(unique) <= 8:
            candidates.append(
                {
                    "column": column,
                    "values": unique[:8],
                    "n_unique": len(unique),
                    "score": round(_score_group_column(column, obs[column]), 2),
                }
            )
    candidates.sort(key=lambda item: item["score"], reverse=True)
    return candidates[:8]


def _score_sample_id_column(column: str, values: pd.Series) -> float:
    text = values.astype("string").fillna("").astype(str).str.strip()
    non_empty = text[text != ""]
    if non_empty.empty:
        return -10.0

    unique = non_empty.unique().tolist()
    unique_count = len(unique)
    row_count = len(non_empty)
    if unique_count < 2:
        return -10.0

    counts = non_empty.value_counts()
    repeated_fraction = float((counts > 1).mean()) if not counts.empty else 0.0
    unique_ratio = unique_count / max(row_count, 1)

    score = 0.0
    canonical = _canonical(column)
    if any(hint in canonical for hint in SAMPLE_ID_COLUMN_HINTS):
        score += 8.0
    if any(hint in canonical for hint in GROUP_COLUMN_HINTS):
        score -= 3.0
    if 2 <= unique_count <= min(96, max(row_count // 4, 2)):
        score += 4.0
    if repeated_fraction >= 0.8:
        score += 3.0
    if counts.min() >= 2:
        score += 1.0
    if unique_ratio >= 0.98:
        score -= 6.0
    return score


def _candidate_sample_id_columns(obs: pd.DataFrame) -> list[dict[str, object]]:
    candidates = []
    for column in obs.columns:
        text = obs[column].astype("string").fillna("").astype(str).str.strip()
        unique = [value for value in text.unique().tolist() if value]
        if len(unique) < 2:
            continue
        score = _score_sample_id_column(column, obs[column])
        if score < 1:
            continue
        candidates.append(
            {
                "column": column,
                "n_unique": len(unique),
                "score": round(score, 2),
            }
        )
    candidates.sort(key=lambda item: item["score"], reverse=True)
    return candidates[:8]


def _choose_pseudobulk_method(values: np.ndarray) -> str:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return "mean"
    if np.nanmin(finite) >= 0 and (np.nanpercentile(finite, 99) > 50 or np.nanmax(finite) > 1000):
        return "sum"
    return "mean"


def _is_likely_single_cell(adata: AnnData, sample_id_candidates: list[dict[str, object]]) -> bool:
    if adata.n_obs < 200:
        return False
    if sample_id_candidates:
        return True
    obs_names = pd.Index(adata.obs_names.astype(str))
    return obs_names.str.contains(r"[-_][ACGT]{6,}", regex=True).mean() >= 0.3


def _split_table_into_matrix_and_metadata(
    indexed: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    numeric = indexed.apply(pd.to_numeric, errors="coerce")
    metadata_columns: list[str] = []

    for column in indexed.columns:
        text = indexed[column].fillna("").astype(str).str.strip()
        unique = [value for value in text.unique().tolist() if value]
        looks_like_metadata = False

        if numeric[column].isna().all():
            looks_like_metadata = True
        elif 2 <= len(unique) <= 8:
            canon = _canonical(column)
            if any(hint in canon for hint in GROUP_COLUMN_HINTS):
                looks_like_metadata = True
            elif _score_group_column(column, indexed[column]) >= 0:
                looks_like_metadata = True

        if looks_like_metadata:
            metadata_columns.append(column)

    metadata = indexed[metadata_columns].copy() if metadata_columns else pd.DataFrame(index=indexed.index)
    matrix = numeric.drop(columns=metadata_columns, errors="ignore")
    matrix = matrix.dropna(axis=0, how="all").dropna(axis=1, how="all")
    return matrix, metadata


def _prepare_h5ad(
    filename: str,
    content: bytes,
    *,
    group_column: str | None = None,
    group1_value: str | None = None,
    group2_value: str | None = None,
    sample_id_column: str | None = None,
) -> PreparedInput:
    with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as temp_file:
        temp_path = Path(temp_file.name)
        temp_file.write(content)
    try:
        adata = load_h5ad_file(str(temp_path))
    finally:
        temp_path.unlink(missing_ok=True)

    values = np.asarray(adata.X)
    if sample_id_column and sample_id_column not in adata.obs.columns:
        raise ValueError(f"Sample identifier column '{sample_id_column}' was not found in adata.obs.")
    candidate_sample_ids = _candidate_sample_id_columns(adata.obs)
    likely_single_cell = _is_likely_single_cell(adata, candidate_sample_ids)
    resolved_sample_id_column = sample_id_column or (
        str(candidate_sample_ids[0]["column"]) if candidate_sample_ids and likely_single_cell else None
    )
    grouping = _auto_group_from_obs(
        adata.obs,
        group_column=group_column,
        group1_value=group1_value,
        group2_value=group2_value,
    )
    if grouping:
        col, value1, value2 = grouping
        if likely_single_cell and resolved_sample_id_column:
            method = _choose_pseudobulk_method(values)
            try:
                bulk = pseudobulk_adata(
                    adata,
                    sample_id_obs_cols=[resolved_sample_id_column],
                    sample_metadata_obs_cols="auto",
                    method=method,
                )
            except Exception as exc:
                raise ValueError(
                    "Single-cell pseudobulk preparation failed objectively. "
                    "Please verify that the selected sample identifier column is consistent within each biological sample."
                ) from exc
            series = compute_vscores_adata(bulk, col, value1, value2)
            prepared = _series_to_adata(series, f"{Path(filename).stem}:{value1}->{value2}")
            status = "ready"
            mode = "single_cell_h5ad_pseudobulk_vscore"
            notes = [
                (
                    f"Detected a single-cell H5AD and pseudobulked cells by '{resolved_sample_id_column}' "
                    f"before computing v-scores from '{col}' using '{value1}' vs '{value2}'."
                ),
                f"Pseudobulk method: {method}.",
            ]
        else:
            series = compute_vscores_adata(adata, col, value1, value2)
            prepared = _series_to_adata(series, f"{Path(filename).stem}:{value1}->{value2}")
            status = "ready"
            mode = "h5ad_auto_vscore"
            notes = [f"Computed v-scores from metadata column '{col}' using '{value1}' vs '{value2}'."]
    elif np.nanmin(values) < 0 and adata.n_obs <= 3:
        prepared = adata.copy()
        status = "ready"
        mode = "signature_h5ad"
        notes = [
            "Detected negative values in a very small H5AD matrix and no objective grouping metadata, so the file was treated as a prepared signature matrix."
        ]
    elif likely_single_cell:
        prepared = None
        status = "needs_configuration"
        mode = "single_cell_h5ad"
        notes = [
            "Detected a likely single-cell H5AD.",
            (
                "An objective DrugReflector signature requires both a two-group column "
                "and a sample identifier column in adata.obs so cells can be pseudobulked first."
            ),
        ]
    else:
        prepared = None
        status = "needs_configuration"
        mode = "h5ad_raw_expression"
        notes = [
            "This H5AD looks like expression data, but no objective two-group contrast could be inferred automatically."
        ]

    summary = _build_summary(
        filename=filename,
        status=status,
        mode=mode,
        prepared=prepared,
        original_shape=[int(adata.n_obs), int(adata.n_vars)],
        notes=notes,
        candidate_groups=_candidate_group_columns(adata.obs),
        candidate_sample_ids=candidate_sample_ids,
        detected_data_type="single_cell_h5ad" if likely_single_cell else "h5ad",
        used_sample_id_column=resolved_sample_id_column,
    )
    return PreparedInput(prepared, summary)


def _prepare_table(
    filename: str,
    content: bytes,
    extension: str,
    *,
    group_column: str | None = None,
    group1_value: str | None = None,
    group2_value: str | None = None,
) -> PreparedInput:
    separator = "\t" if extension == ".tsv" else ","
    table = pd.read_csv(io.BytesIO(content), sep=separator)
    if table.empty:
        raise ValueError("Uploaded table is empty.")

    # Differential result table: gene + score/stat column.
    gene_column = _choose_gene_column(table.columns.tolist())
    score_column = _choose_score_column(table)
    gene_ratio = (
        _gene_likeness(pd.Index(table[gene_column].dropna().astype(str)))
        if gene_column
        else 0.0
    )
    if gene_column and score_column and gene_ratio >= 0.5 and len(table.columns) <= 12:
        series = pd.Series(table[score_column].to_numpy(dtype=float), index=table[gene_column])
        prepared = _series_to_adata(series, Path(filename).stem or "signature")
        summary = _build_summary(
            filename=filename,
            status="ready",
            mode="differential_table",
            prepared=prepared,
            original_shape=[int(table.shape[0]), int(table.shape[1])],
            notes=[f"Converted differential table using gene column '{gene_column}' and score column '{score_column}'."],
        )
        return PreparedInput(prepared, summary)

    indexed = table.set_index(table.columns[0])
    matrix, sample_metadata = _split_table_into_matrix_and_metadata(indexed)
    if matrix.empty:
        raise ValueError("Could not find a numeric matrix in the uploaded table.")

    row_gene_score = _gene_likeness(matrix.index)
    col_gene_score = _gene_likeness(matrix.columns)

    if row_gene_score > col_gene_score:
        sample_by_gene = matrix.T
        orientation = "gene_by_sample"
    else:
        sample_by_gene = matrix
        orientation = "sample_by_gene"

    sample_by_gene = sample_by_gene.fillna(0.0)
    values = sample_by_gene.to_numpy(dtype=float, copy=False)

    if orientation == "sample_by_gene" and not sample_metadata.empty and sample_by_gene.shape[0] > 3:
        sample_metadata = sample_metadata.reindex(sample_by_gene.index)
        grouping = _auto_group_from_obs(
            sample_metadata,
            group_column=group_column,
            group1_value=group1_value,
            group2_value=group2_value,
        )
        if grouping:
            col, value1, value2 = grouping
            adata = AnnData(
                X=sample_by_gene.to_numpy(dtype=np.float32, copy=True),
                obs=sample_metadata.copy(),
                var=pd.DataFrame(index=[_normalize_gene_name(col) for col in sample_by_gene.columns]),
            )
            series = compute_vscores_adata(adata, col, value1, value2)
            prepared = _series_to_adata(series, f"{Path(filename).stem}:{value1}->{value2}")
            summary = _build_summary(
                filename=filename,
                status="ready",
                mode="expression_matrix_metadata_vscore",
                prepared=prepared,
                original_shape=[int(indexed.shape[0]), int(indexed.shape[1])],
                notes=[f"Computed v-scores from metadata column '{col}' using '{value1}' vs '{value2}'."],
                candidate_groups=_candidate_group_columns(sample_metadata),
            )
            return PreparedInput(prepared, summary)

    if np.nanmin(values) < 0 and sample_by_gene.shape[0] <= 3:
        prepared = _frame_to_adata(sample_by_gene)
        summary = _build_summary(
            filename=filename,
            status="ready",
            mode="signature_matrix",
            prepared=prepared,
            original_shape=[int(indexed.shape[0]), int(indexed.shape[1])],
            notes=[f"Detected a prepared signature matrix ({orientation.replace('_', ' ')})."],
            candidate_groups=_candidate_group_columns(sample_metadata) if orientation == "sample_by_gene" else None,
        )
        return PreparedInput(prepared, summary)

    label_groups = _build_group_series_from_labels(sample_by_gene.index)
    if label_groups is not None and _is_raw_expression(values):
        adata = AnnData(
            X=sample_by_gene.to_numpy(dtype=np.float32, copy=True),
            obs=pd.DataFrame({"_label_group": label_groups}, index=sample_by_gene.index.astype(str)),
            var=pd.DataFrame(index=[_normalize_gene_name(col) for col in sample_by_gene.columns]),
        )
        groups = list(pd.Index(label_groups).unique())
        series = compute_vscores_adata(adata, "_label_group", groups[0], groups[1])
        prepared = _series_to_adata(series, f"{Path(filename).stem}:{groups[0]}->{groups[1]}")
        summary = _build_summary(
            filename=filename,
            status="ready",
            mode="expression_matrix_auto_vscore",
            prepared=prepared,
            original_shape=[int(indexed.shape[0]), int(indexed.shape[1])],
            notes=[f"Computed v-scores from sample labels using inferred groups '{groups[0]}' vs '{groups[1]}'."],
        )
        return PreparedInput(prepared, summary)

    if _is_raw_expression(values) and sample_by_gene.shape[0] > 3:
        summary = _build_summary(
            filename=filename,
            status="needs_configuration",
            mode="raw_expression_matrix",
            prepared=None,
            original_shape=[int(indexed.shape[0]), int(indexed.shape[1])],
            notes=[
                f"Interpreted the table as a {orientation.replace('_', ' ')} matrix.",
                "Values look like expression data, but no objective two-group contrast could be inferred from metadata or sample labels.",
            ],
            candidate_groups=_candidate_group_columns(sample_metadata) if orientation == "sample_by_gene" else None,
        )
        return PreparedInput(None, summary)

    if np.nanmin(values) < 0:
        summary = _build_summary(
            filename=filename,
            status="needs_configuration",
            mode="ambiguous_negative_matrix",
            prepared=None,
            original_shape=[int(indexed.shape[0]), int(indexed.shape[1])],
            notes=[
                f"Interpreted the table as a {orientation.replace('_', ' ')} matrix.",
                "Negative values were detected, but the file contains too many samples to treat as a prepared signature matrix without an objective contrast definition.",
            ],
            candidate_groups=_candidate_group_columns(sample_metadata) if orientation == "sample_by_gene" else None,
        )
        return PreparedInput(None, summary)

    if sample_by_gene.shape[0] <= 3:
        prepared = _frame_to_adata(sample_by_gene)
        summary = _build_summary(
            filename=filename,
            status="ready",
            mode="signature_matrix",
            prepared=prepared,
            original_shape=[int(indexed.shape[0]), int(indexed.shape[1])],
            notes=[f"Interpreted the table as a {orientation.replace('_', ' ')} matrix with a small number of samples."],
            candidate_groups=_candidate_group_columns(sample_metadata) if orientation == "sample_by_gene" else None,
        )
        return PreparedInput(prepared, summary)

    summary = _build_summary(
        filename=filename,
        status="needs_configuration",
        mode="unresolved_matrix",
        prepared=None,
        original_shape=[int(indexed.shape[0]), int(indexed.shape[1])],
        notes=[
            f"Interpreted the table as a {orientation.replace('_', ' ')} matrix.",
            "The file could not be converted into an objective signature automatically.",
        ],
        candidate_groups=_candidate_group_columns(sample_metadata) if orientation == "sample_by_gene" else None,
    )
    return PreparedInput(None, summary)


def _build_summary(
    *,
    filename: str,
    status: str,
    mode: str,
    prepared: AnnData | None,
    original_shape: list[int],
    notes: list[str],
    candidate_groups: list[dict[str, object]] | None = None,
    candidate_sample_ids: list[dict[str, object]] | None = None,
    detected_data_type: str | None = None,
    used_sample_id_column: str | None = None,
) -> dict[str, object]:
    summary: dict[str, object] = {
        "filename": filename,
        "status": status,
        "mode": mode,
        "original_shape": original_shape,
        "notes": notes,
        "candidate_groups": candidate_groups or [],
        "candidate_sample_ids": candidate_sample_ids or [],
        "detected_data_type": detected_data_type,
        "used_sample_id_column": used_sample_id_column,
    }
    if prepared is None:
        summary["prepared_shape"] = None
        summary["sample_names"] = []
        summary["top_genes"] = []
        return summary

    scores = pd.Series(prepared.X[0], index=prepared.var_names) if prepared.n_obs == 1 else None
    top_genes = []
    if scores is not None:
        order = scores.abs().sort_values(ascending=False).head(8).index
        top_genes = [
            {"gene": gene, "score": float(scores.loc[gene])}
            for gene in order
        ]

    summary["prepared_shape"] = [int(prepared.n_obs), int(prepared.n_vars)]
    summary["sample_names"] = prepared.obs_names.tolist()[:8]
    summary["top_genes"] = top_genes
    summary["n_negative_values"] = int(np.sum(np.asarray(prepared.X) < 0))
    return summary


def prepare_uploaded_input(
    filename: str,
    content: bytes,
    *,
    group_column: str | None = None,
    group1_value: str | None = None,
    group2_value: str | None = None,
    sample_id_column: str | None = None,
) -> PreparedInput:
    extension = Path(filename).suffix.lower()
    if extension == ".h5ad":
        return _prepare_h5ad(
            filename,
            content,
            group_column=group_column,
            group1_value=group1_value,
            group2_value=group2_value,
            sample_id_column=sample_id_column,
        )
    if extension in {".csv", ".tsv"}:
        return _prepare_table(
            filename,
            content,
            extension,
            group_column=group_column,
            group1_value=group1_value,
            group2_value=group2_value,
        )
    raise ValueError(f"Unsupported file type: {extension!r}")
