"""Objective direction-vs-signature evidence helpers."""

from __future__ import annotations

from functools import lru_cache
import json
import os
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from anndata import AnnData

CONNECTIVITY_TABLE_ENV_KEYS = (
    "DRUGREFLECTOR_CONNECTIVITY_TABLE",
    "CLUE_CONNECTIVITY_TABLE",
    "LINCS_CONNECTIVITY_TABLE",
)
DEFAULT_SOURCE = "External signed connectivity table"
L1000CDS2_URL = "https://maayanlab.cloud/L1000CDS2/query"


def _to_signature_frame(data: pd.DataFrame | AnnData) -> pd.DataFrame:
    if isinstance(data, pd.DataFrame):
        return data.copy()

    return pd.DataFrame(
        np.asarray(data.X, dtype=float),
        index=data.obs_names.astype(str),
        columns=data.var_names.astype(str),
    )


def _core_brd_id(value: object) -> str:
    text = str(value or "").strip().upper()
    match = pd.Series([text]).str.extract(r"((?:BRD|BRDN)-[A-Z0-9]+)", expand=False).iloc[0]
    return str(match) if pd.notna(match) else text


def _normalize_alias(value: object) -> str:
    text = str(value or "").strip().lower()
    return "".join(ch for ch in text if ch.isalnum())


def _normalize_column_name(value: object) -> str:
    return "".join(ch for ch in str(value).strip().lower() if ch.isalnum())


def _first_matching_column(columns: pd.Index, candidates: set[str]) -> str | None:
    normalized = {_normalize_column_name(column): str(column) for column in columns}
    for candidate in candidates:
        if candidate in normalized:
            return normalized[candidate]
    return None


def _connectivity_table_path() -> Path | None:
    for key in CONNECTIVITY_TABLE_ENV_KEYS:
        value = os.getenv(key)
        if value:
            return Path(value).expanduser()
    return None


def _read_table(path: Path) -> pd.DataFrame:
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".tsv") or suffixes.endswith(".txt"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


@lru_cache(maxsize=1)
def _load_connectivity_table() -> tuple[pd.DataFrame | None, str | None]:
    path = _connectivity_table_path()
    if path is None:
        return None, None
    if not path.exists():
        return None, f"Configured signed connectivity table does not exist: {path}"

    try:
        raw = _read_table(path)
    except Exception as exc:
        return None, f"Could not read signed connectivity table '{path}': {exc}"

    compound_col = _first_matching_column(
        raw.columns,
        {
            "compound",
            "pertid",
            "brdid",
            "broadid",
            "perturbagen",
            "perturbagenid",
        },
    )
    score_col = _first_matching_column(
        raw.columns,
        {
            "connectivityscore",
            "signedconnectivityscore",
            "score",
            "cs",
            "tau",
            "normconnectivityscore",
            "normalizedconnectivityscore",
        },
    )
    if not compound_col or not score_col:
        return (
            None,
            "Signed connectivity table must contain a compound/pert_id/BRD column "
            "and a signed score column such as connectivity_score, tau, cs, or score.",
        )

    sample_col = _first_matching_column(
        raw.columns,
        {"signature", "sample", "query", "queryname", "inputsignature"},
    )
    source_col = _first_matching_column(raw.columns, {"source", "dataset", "evidence"})

    table = pd.DataFrame(
        {
            "compound": raw[compound_col].map(_core_brd_id),
            "score": pd.to_numeric(raw[score_col], errors="coerce"),
            "sample": raw[sample_col].astype(str) if sample_col else "",
            "source": raw[source_col].astype(str) if source_col else DEFAULT_SOURCE,
        }
    )
    table = table[table["compound"].ne("") & table["score"].notna()].copy()
    if table.empty:
        return None, f"Signed connectivity table '{path}' did not contain usable signed scores."

    return table, None


def _score_threshold(table: pd.DataFrame | None) -> float:
    configured = os.getenv("DRUGREFLECTOR_CONNECTIVITY_SCORE_THRESHOLD")
    if configured:
        try:
            return abs(float(configured))
        except ValueError:
            pass
    if table is not None and not table.empty and table["score"].abs().max() <= 1.5:
        return 0.1
    return 20.0


def _direction_from_score(score: float, threshold: float) -> str:
    if score <= -threshold:
        return "Reverse"
    if score >= threshold:
        return "Mimic"
    return "No objective evidence"


def _query_gene_counts(values: pd.Series) -> tuple[int, int]:
    return int((values > 0).sum()), int((values < 0).sum())


def _fallback_context(values: pd.Series, reason: str, *, source: str | None = None) -> dict[str, Any]:
    up_count, down_count = _query_gene_counts(values)
    return {
        "label": "No objective evidence",
        "score": None,
        "source": source or "Signed connectivity evidence not available",
        "reason": reason,
        "query_up_genes": up_count,
        "query_down_genes": down_count,
    }


def _external_l1000cds2_enabled() -> bool:
    value = os.getenv("DRUGREFLECTOR_L1000CDS2_ENABLED", "1").strip().lower()
    return value not in {"0", "false", "no", "off"}


def _l1000cds2_timeout() -> float:
    value = os.getenv("DRUGREFLECTOR_L1000CDS2_TIMEOUT_SECONDS", "20").strip()
    try:
        return max(3.0, float(value))
    except ValueError:
        return 20.0


def _l1000_query_gene_count() -> int:
    value = os.getenv("DRUGREFLECTOR_L1000CDS2_TOP_GENES", "150").strip()
    try:
        return max(10, min(500, int(value)))
    except ValueError:
        return 150


def _l1000_score_threshold() -> float:
    value = os.getenv("DRUGREFLECTOR_L1000CDS2_SCORE_THRESHOLD", "0").strip()
    try:
        return max(0.0, float(value))
    except ValueError:
        return 0.0


def _signature_query_genes(values: pd.Series) -> tuple[tuple[str, ...], tuple[str, ...]]:
    limit = _l1000_query_gene_count()
    ordered = values.dropna().sort_values(ascending=False)
    up = tuple(ordered[ordered > 0].head(limit).index.astype(str).str.upper())
    down = tuple(ordered[ordered < 0].sort_values(ascending=True).head(limit).index.astype(str).str.upper())
    return up, down


@lru_cache(maxsize=256)
def _query_l1000cds2(up_genes: tuple[str, ...], down_genes: tuple[str, ...], *, mimic: bool) -> tuple[dict[str, Any] | None, str | None]:
    if not _external_l1000cds2_enabled():
        return None, "L1000CDS2 querying is disabled by DRUGREFLECTOR_L1000CDS2_ENABLED."
    if not up_genes or not down_genes:
        return None, "L1000CDS2 query requires both up and down gene sets."

    payload = {
        "data": {"upGenes": list(up_genes), "dnGenes": list(down_genes)},
        "config": {
            "aggravate": mimic,
            "searchMethod": "geneSet",
            "share": False,
            "combination": False,
            "db-version": os.getenv("DRUGREFLECTOR_L1000CDS2_DB_VERSION", "latest"),
        },
        "meta": [{"key": "Tag", "value": "DrugReflector signed direction query"}],
    }

    try:
        request = Request(
            L1000CDS2_URL,
            data=json.dumps(payload).encode("utf-8"),
            headers={"content-type": "application/json"},
            method="POST",
        )
        with urlopen(request, timeout=_l1000cds2_timeout()) as response:
            raw = response.read().decode("utf-8", errors="replace")
        return json.loads(raw), None
    except HTTPError as exc:
        return None, f"L1000CDS2 query failed with HTTP {exc.code}: {exc.reason}"
    except URLError as exc:
        return None, f"L1000CDS2 query failed: {exc.reason}"
    except Exception as exc:
        return None, f"L1000CDS2 query failed: {exc}"


def _compound_aliases(compound: str, aliases: dict[str, list[str]] | None) -> set[str]:
    values = {compound, _core_brd_id(compound)}
    if aliases:
        values.update(aliases.get(str(compound), []))
        values.update(aliases.get(_core_brd_id(compound), []))
    normalized = {_normalize_alias(value) for value in values if str(value).strip()}
    normalized.discard("")
    return normalized


def _match_l1000_results(result: dict[str, Any] | None, compound: str, aliases: dict[str, list[str]] | None) -> dict[str, Any] | None:
    if not result:
        return None
    targets = _compound_aliases(compound, aliases)
    rows = result.get("topMeta") or result.get("results") or []
    if not isinstance(rows, list):
        return None

    matches: list[dict[str, Any]] = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        candidates = {
            _core_brd_id(row.get("pert_id")),
            str(row.get("pert_desc") or ""),
            str(row.get("pert_iname") or ""),
            str(row.get("compound") or ""),
        }
        normalized = {_normalize_alias(value) for value in candidates if str(value).strip()}
        if targets & normalized:
            matches.append(row)

    if not matches:
        return None
    return max(matches, key=lambda item: float(item.get("score") or 0.0))


def _l1000_direction_context(
    values: pd.Series,
    *,
    sample_name: str,
    compound: str,
    aliases: dict[str, list[str]] | None,
) -> dict[str, Any] | None:
    up_genes, down_genes = _signature_query_genes(values)
    reverse_result, reverse_error = _query_l1000cds2(up_genes, down_genes, mimic=False)
    mimic_result, mimic_error = _query_l1000cds2(up_genes, down_genes, mimic=True)

    reverse_match = _match_l1000_results(reverse_result, compound, aliases)
    mimic_match = _match_l1000_results(mimic_result, compound, aliases)
    threshold = _l1000_score_threshold()

    candidates: list[tuple[str, dict[str, Any]]] = []
    if reverse_match and float(reverse_match.get("score") or 0.0) >= threshold:
        candidates.append(("Reverse", reverse_match))
    if mimic_match and float(mimic_match.get("score") or 0.0) >= threshold:
        candidates.append(("Mimic", mimic_match))

    if candidates:
        label, row = max(candidates, key=lambda item: float(item[1].get("score") or 0.0))
        score = float(row.get("score") or 0.0)
        description = row.get("pert_desc") or row.get("pert_id") or compound
        cell = row.get("cell_id")
        dose = row.get("pert_dose")
        time = row.get("pert_time")
        reason = (
            f"Matched real L1000CDS2 {label.lower()} result for {description} "
            f"(score={score:.4g}, cell={cell}, dose={dose}, time={time})."
        )
        up_count, down_count = _query_gene_counts(values)
        return {
            "label": label,
            "score": score if label == "Mimic" else -score,
            "source": "L1000CDS2 LINCS L1000 gene-set query",
            "reason": reason,
            "query_up_genes": up_count,
            "query_down_genes": down_count,
            "external_signature_id": row.get("sig_id"),
            "external_pert_id": row.get("pert_id"),
            "external_pert_desc": row.get("pert_desc"),
            "external_cell_id": row.get("cell_id"),
        }

    if reverse_error and mimic_error:
        return _fallback_context(
            values,
            f"No objective evidence because real L1000CDS2 reverse and mimic queries failed. "
            f"Reverse error: {reverse_error}; Mimic error: {mimic_error}",
            source="L1000CDS2 LINCS L1000 gene-set query",
        )

    return _fallback_context(
        values,
        "No objective evidence because real L1000CDS2 reverse/mimic queries completed, "
        "but this DrugReflector compound was not found among the returned top LINCS perturbation signatures.",
        source="L1000CDS2 LINCS L1000 gene-set query",
    )


def _sample_direction_context(values: pd.Series, *, min_genes: int = 10) -> dict[str, Any]:
    ordered = values.sort_values(ascending=False)
    up = ordered[ordered > 0]
    down = ordered[ordered < 0]

    if len(up) < min_genes or len(down) < min_genes:
        return _fallback_context(
            values,
            "No objective evidence because the prepared input signature does not contain at least "
            f"{min_genes} positive and {min_genes} negative genes for a bidirectional connectivity query.",
        )

    table, error = _load_connectivity_table()
    if error:
        return _fallback_context(values, error)
    if table is None:
        return _fallback_context(
            values,
            "No objective evidence because no real CLUE/LINCS signed connectivity result table is configured. "
            "Set DRUGREFLECTOR_CONNECTIVITY_TABLE to a CSV/TSV exported from CLUE/LINCS before classifying "
            "compounds as Reverse or Mimic.",
        )

    return _fallback_context(
        values,
        "No objective evidence because the configured signed connectivity table does not contain this compound/signature pair.",
        source=DEFAULT_SOURCE,
    )


def _compound_direction_context(
    values: pd.Series,
    *,
    sample_name: str,
    compound: str,
    table: pd.DataFrame | None,
    table_error: str | None,
    compound_aliases: dict[str, list[str]] | None = None,
    min_genes: int = 10,
) -> dict[str, Any]:
    up_count, down_count = _query_gene_counts(values)
    if up_count < min_genes or down_count < min_genes:
        return _sample_direction_context(values, min_genes=min_genes)

    if table_error:
        return _fallback_context(values, table_error)
    if table is None:
        l1000_context = _l1000_direction_context(
            values,
            sample_name=sample_name,
            compound=compound,
            aliases=compound_aliases,
        )
        if l1000_context is not None:
            return l1000_context
        return _sample_direction_context(values, min_genes=min_genes)

    compound_key = _core_brd_id(compound)
    subset = table[table["compound"].eq(compound_key)].copy()
    if subset.empty:
        return _fallback_context(
            values,
            "No objective evidence because this compound was not found in the configured real signed connectivity table.",
            source=DEFAULT_SOURCE,
        )

    exact = subset[subset["sample"].eq(str(sample_name))]
    if not exact.empty:
        subset = exact
    elif subset["sample"].ne("").any():
        global_subset = subset[subset["sample"].eq("")]
        if not global_subset.empty:
            subset = global_subset
        else:
            return _fallback_context(
                values,
                "No objective evidence because this compound exists in the signed connectivity table, "
                "but not for the current input signature/query name.",
                source=DEFAULT_SOURCE,
            )

    row = subset.iloc[subset["score"].abs().argmax()]
    score = float(row["score"])
    threshold = _score_threshold(table)
    label = _direction_from_score(score, threshold)
    if label == "Reverse":
        reason = (
            f"Real signed connectivity score {score:.4g} is <= -{threshold:g}; "
            "the compound perturbation is objectively opposite to the input signature."
        )
    elif label == "Mimic":
        reason = (
            f"Real signed connectivity score {score:.4g} is >= {threshold:g}; "
            "the compound perturbation is objectively concordant with the input signature."
        )
    else:
        reason = (
            f"Real signed connectivity score {score:.4g} is below the objective threshold "
            f"(|score| < {threshold:g}), so the direction is not classified."
        )

    return {
        "label": label,
        "score": score,
        "source": str(row.get("source") or DEFAULT_SOURCE),
        "reason": reason,
        "query_up_genes": up_count,
        "query_down_genes": down_count,
    }


def get_direction_evidence(
    data: pd.DataFrame | AnnData,
    compounds_by_sample: dict[str, list[str]],
    compound_aliases: dict[str, list[str]] | None = None,
) -> dict[str, dict[str, dict[str, Any]]]:
    frame = _to_signature_frame(data)
    evidence: dict[str, dict[str, dict[str, Any]]] = {}
    table, table_error = _load_connectivity_table()

    for sample_name, compounds in compounds_by_sample.items():
        if sample_name not in frame.index:
            continue
        evidence[sample_name] = {
            str(compound): _compound_direction_context(
                frame.loc[sample_name],
                sample_name=str(sample_name),
                compound=str(compound),
                table=table,
                table_error=table_error,
                compound_aliases=compound_aliases,
            )
            for compound in compounds
        }

    return evidence
