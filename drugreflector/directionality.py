"""Objective direction-vs-signature evidence helpers."""

from __future__ import annotations

import os
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData


def _to_signature_frame(data: pd.DataFrame | AnnData) -> pd.DataFrame:
    if isinstance(data, pd.DataFrame):
        return data.copy()

    return pd.DataFrame(
        np.asarray(data.X, dtype=float),
        index=data.obs_names.astype(str),
        columns=data.var_names.astype(str),
    )


def _sample_direction_context(values: pd.Series, *, min_genes: int = 10) -> dict[str, Any]:
    ordered = values.sort_values(ascending=False)
    up = ordered[ordered > 0]
    down = ordered[ordered < 0]

    if len(up) < min_genes or len(down) < min_genes:
        reason = (
            "No objective evidence because the prepared input signature does not contain at least "
            f"{min_genes} positive and {min_genes} negative genes for a bidirectional connectivity query."
        )
    elif not os.getenv("CLUE_API_KEY"):
        reason = (
            "No objective evidence because this deployment has no configured CLUE Touchstone API key, "
            "so an external signed connectivity query cannot be run objectively."
        )
    else:
        reason = (
            "No objective evidence because automated CLUE Touchstone querying is not enabled in this build."
        )

    return {
        "label": "No objective evidence",
        "score": None,
        "source": "CLUE Touchstone connectivity (not available in this deployment)",
        "reason": reason,
        "query_up_genes": int(len(up)),
        "query_down_genes": int(len(down)),
    }


def get_direction_evidence(
    data: pd.DataFrame | AnnData,
    compounds_by_sample: dict[str, list[str]],
) -> dict[str, dict[str, dict[str, Any]]]:
    frame = _to_signature_frame(data)
    evidence: dict[str, dict[str, dict[str, Any]]] = {}

    for sample_name, compounds in compounds_by_sample.items():
        if sample_name not in frame.index:
            continue
        sample_context = _sample_direction_context(frame.loc[sample_name])
        evidence[sample_name] = {
            str(compound): dict(sample_context)
            for compound in compounds
        }

    return evidence
