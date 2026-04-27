"""
Utility functions for DrugReflector.

This module contains the preprocessing helpers that are expected to be
available from the installed package, including v-score computation,
AnnData loading, and pseudobulking utilities.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from anndata import AnnData, read_h5ad
from scipy.optimize import minimize_scalar


def _norm_func(data, clip, target_std):
    """Helper function called in `clip_rescale_rows`."""
    return lambda norm: (np.std(np.clip(data / norm, -clip, clip)) - target_std) ** 2


def clip_rescale_rows(X, clip, target_std, bounds=(0, 1e3)):
    """Rescale and clip each row of ``X`` to match the target standard deviation."""
    n = X.shape[0]
    norm = np.zeros(n)
    for i in range(n):
        ms = minimize_scalar(
            _norm_func(X[i], clip, target_std),
            method="Bounded",
            bounds=bounds,
        )
        assert ms.success, f"Unable to rescale row {i}."
        norm[i] = ms.x

    np.divide(X, norm[:, np.newaxis], out=X)
    np.clip(X, -clip, clip, out=X)


def compute_vscores(adata, transitions=None, mute=False):
    """
    Compute v-scores for gene expression data.

    If ``transitions`` is omitted, the input AnnData is assumed to already
    contain v-scores and is returned unchanged.
    """
    if transitions is not None:
        if not isinstance(transitions, dict):
            raise ValueError(
                "transitions must be a dict with group_col, group1_value, group2_value keys"
            )

        required_keys = ["group_col", "group1_value", "group2_value"]
        for key in required_keys:
            if key not in transitions:
                raise ValueError(f"transitions dict must contain key: {key}")

        group_col = transitions["group_col"]
        group1_value = transitions["group1_value"]
        group2_value = transitions["group2_value"]
        layer = transitions.get("layer", None)

        vscores_series = compute_vscores_adata(
            adata,
            group_col,
            group1_value,
            group2_value,
            layer=layer,
        )

        obs_name = vscores_series.name
        return AnnData(
            X=vscores_series.values.reshape(1, -1),
            var=pd.DataFrame(index=vscores_series.index),
            obs=pd.DataFrame(index=[obs_name]),
        )

    if not mute:
        warnings.warn("Assuming passed representation is v-score.", stacklevel=1)
    return adata.copy()


def load_h5ad_file(filepath: str) -> AnnData:
    """Load an H5AD file and sanitize its matrix for inference."""
    adata = read_h5ad(filepath)

    if adata.X is None:
        raise ValueError("AnnData object has no X matrix")

    if hasattr(adata.X, "toarray"):
        adata.X = adata.X.toarray()

    if not np.isfinite(adata.X).all():
        print("Warning: Found infinite or NaN values in data, replacing with 0")
        adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

    return adata


def create_synthetic_gene_expression(
    n_obs: int,
    n_vars: int,
    obs_names: list | None = None,
    var_names: list | None = None,
    random_state: int = 42,
) -> AnnData:
    """Create synthetic expression data for local testing."""
    np.random.seed(random_state)
    X = np.random.normal(0, 1, size=(n_obs, n_vars))

    if obs_names is None:
        obs_names = [f"sample_{i}" for i in range(n_obs)]

    if var_names is None:
        var_names = [f"gene_{i}" for i in range(n_vars)]

    return AnnData(
        X=X,
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=var_names),
    )


def pseudobulk_adata(
    adata,
    sample_id_obs_cols,
    sample_metadata_obs_cols="auto",
    layer=None,
    method="sum",
):
    """
    Pseudobulk an AnnData object by aggregating rows over sample identifiers.
    """
    import scipy.sparse as sp

    adata_temp = adata.copy()

    if layer:
        if layer not in adata_temp.layers:
            raise ValueError(f"Layer '{layer}' not found in AnnData object")
        X_data = adata_temp.layers[layer]
    else:
        X_data = adata_temp.X

    is_sparse = sp.issparse(X_data)
    if is_sparse:
        X_data = X_data.tocsr()

    adata_obs_cols = set(adata_temp.obs.columns)
    adata_obs_cols -= set(sample_id_obs_cols)

    adata_temp.obs["_TempIndex"] = adata_temp.obs.apply(
        lambda row: "_".join([f"{row[id_col]}" for id_col in sample_id_obs_cols]),
        axis="columns",
    )

    group_names = adata_temp.obs["_TempIndex"].unique()
    n_groups = len(group_names)
    n_genes = adata_temp.n_vars

    if is_sparse:
        bulk_X = sp.lil_matrix((n_groups, n_genes), dtype=X_data.dtype)
    else:
        bulk_X = np.zeros((n_groups, n_genes), dtype=X_data.dtype)

    for i, group_name in enumerate(group_names):
        group_mask = adata_temp.obs["_TempIndex"] == group_name
        group_data = X_data[group_mask]

        if method == "sum":
            bulk_X[i, :] = group_data.sum(axis=0)
        elif method == "mean":
            bulk_X[i, :] = group_data.mean(axis=0)
        else:
            raise ValueError("method parameter must be 'sum' or 'mean'")

    if is_sparse:
        bulk_X = bulk_X.tocsr()

    bulk_adata = AnnData(
        X=bulk_X,
        var=adata_temp.var.copy(),
        dtype=X_data.dtype,
    )
    bulk_adata.obs.index = group_names

    if sample_metadata_obs_cols == "auto":
        sample_metadata_obs_cols = []
        for obs_col in adata_obs_cols:
            if adata_temp.obs.groupby("_TempIndex")[obs_col].nunique().max() == 1:
                sample_metadata_obs_cols.append(obs_col)

    if sample_metadata_obs_cols or sample_id_obs_cols:
        try:
            cols_to_use = ["_TempIndex"] + list(sample_id_obs_cols)
            if sample_metadata_obs_cols:
                cols_to_use.extend(sample_metadata_obs_cols)

            metadata_mapping = (
                adata_temp.obs[cols_to_use]
                .drop_duplicates()
                .set_index("_TempIndex", verify_integrity=True)
            )
        except ValueError as exc:
            raise ValueError(
                "The combination of values in sample_metadata_cols of adata.obs "
                "must be unique for each value in sample_id_col."
            ) from exc

        bulk_adata.obs = bulk_adata.obs.merge(
            metadata_mapping,
            left_index=True,
            right_index=True,
            how="left",
        )

    if hasattr(adata_temp, "layers") and adata_temp.layers:
        for layer_name, layer_data in adata_temp.layers.items():
            if layer_name == layer:
                continue

            is_layer_sparse = sp.issparse(layer_data)
            if is_layer_sparse:
                layer_data = layer_data.tocsr()

            if is_layer_sparse:
                layer_bulk = sp.lil_matrix((n_groups, n_genes), dtype=layer_data.dtype)
            else:
                layer_bulk = np.zeros((n_groups, n_genes), dtype=layer_data.dtype)

            for i, group_name in enumerate(group_names):
                group_mask = adata_temp.obs["_TempIndex"] == group_name
                group_layer_data = layer_data[group_mask]

                if method == "sum":
                    layer_bulk[i, :] = group_layer_data.sum(axis=0)
                elif method == "mean":
                    layer_bulk[i, :] = group_layer_data.mean(axis=0)

            if is_layer_sparse:
                layer_bulk = layer_bulk.tocsr()

            bulk_adata.layers[layer_name] = layer_bulk

    return bulk_adata


def compute_vscore_two_groups(group1, group2):
    """Compute the v-score between two groups of values."""
    group1 = np.asarray(group1)
    group2 = np.asarray(group2)

    mean1 = np.mean(group1)
    mean2 = np.mean(group2)
    var1 = np.var(group1, ddof=0)
    var2 = np.var(group2, ddof=0)

    denominator = np.sqrt(var1 + var2) + (var1 + var2 == 0)
    return (mean2 - mean1) / denominator


def compute_vscores_adata(adata, group_col, group1_value, group2_value, layer=None):
    """Compute vectorized v-scores between two populations in an AnnData object."""
    import scipy.sparse as sp

    if group_col not in adata.obs.columns:
        raise ValueError(f"Column '{group_col}' not found in adata.obs")

    group1_mask = adata.obs[group_col] == group1_value
    group2_mask = adata.obs[group_col] == group2_value

    if not group1_mask.any():
        raise ValueError(f"No samples found with {group_col}='{group1_value}'")
    if not group2_mask.any():
        raise ValueError(f"No samples found with {group_col}='{group2_value}'")

    if layer and layer in adata.layers:
        X_data = adata.layers[layer]
    else:
        X_data = adata.X

    if sp.issparse(X_data):
        X_data = X_data.toarray()

    group1_data = X_data[group1_mask, :]
    group2_data = X_data[group2_mask, :]

    mean1 = np.mean(group1_data, axis=0)
    mean2 = np.mean(group2_data, axis=0)
    var1 = np.var(group1_data, axis=0, ddof=0)
    var2 = np.var(group2_data, axis=0, ddof=0)

    denominator = np.sqrt(var1 + var2) + (var1 + var2 == 0)
    vscores = (mean2 - mean1) / denominator

    series_name = f"{group_col}:{group1_value}->{group2_value}"
    return pd.Series(vscores, index=adata.var_names, name=series_name)
