"""
Package-local signature refinement implementation.

Keeping this module inside the installed package avoids import ambiguity when
the project is used through editable installs or console entrypoints.
"""

from __future__ import annotations

import warnings
from typing import List, Optional, Union

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData, concat
from scipy.stats import pearsonr

from .utils import pseudobulk_adata


class SignatureRefinement:
    """
    Refine transcriptional signatures using paired transcriptional + phenotypic data.
    """

    def __init__(self, starting_signature: Union[AnnData, pd.Series]):
        if isinstance(starting_signature, AnnData):
            if starting_signature.n_obs != 1:
                raise ValueError(
                    "AnnData starting signature must have exactly 1 observation (row)"
                )
            self.starting_signature = pd.Series(
                starting_signature.X.flatten(),
                index=starting_signature.var_names,
                name="starting_signature",
            )
        elif isinstance(starting_signature, pd.Series):
            self.starting_signature = starting_signature.copy()
        else:
            raise ValueError("starting_signature must be AnnData or pd.Series")

        self.expr = None
        self.readouts = None
        self.learned_signatures = None
        self.refined_signatures = None
        self._compound_id_obs_col = None
        self._sample_id_obs_cols = None
        self._signature_id_obs_cols = None

    def load_counts_data(
        self,
        adata: AnnData,
        compound_id_obs_col: str,
        layer: Optional[str] = None,
        sample_id_obs_cols: Optional[List[str]] = None,
        signature_id_obs_cols: Optional[List[str]] = None,
    ):
        if layer and layer in adata.layers:
            counts_data = adata.layers[layer]
        else:
            counts_data = adata.X

        if hasattr(counts_data, "data"):
            sample_values = counts_data.data[:1000]
        else:
            sample_values = counts_data.flatten()[:1000]

        if not np.all(sample_values >= 0):
            warnings.warn("Data contains negative values - may not be raw counts!")
        elif not np.allclose(sample_values, sample_values.astype(int)):
            warnings.warn("Data contains non-integer values - may not be raw counts!")

        if sample_id_obs_cols is None:
            sample_id_obs_cols = []
        if signature_id_obs_cols is None:
            signature_id_obs_cols = []

        sample_pseudobulk_cols = (
            [compound_id_obs_col] + sample_id_obs_cols + signature_id_obs_cols
        )
        pb_adata = pseudobulk_adata(
            adata,
            sample_id_obs_cols=sample_pseudobulk_cols,
            layer=layer,
            method="sum",
        )

        if not hasattr(pb_adata, "layers"):
            pb_adata.layers = {}
        pb_adata.layers["pseudobulked_counts"] = pb_adata.X.copy()

        sc.pp.normalize_total(pb_adata, target_sum=1e6)
        sc.pp.log1p(pb_adata)

        self.expr = pb_adata
        self._compound_id_obs_col = compound_id_obs_col
        self._sample_id_obs_cols = sample_id_obs_cols if sample_id_obs_cols else []
        self._signature_id_obs_cols = (
            signature_id_obs_cols if signature_id_obs_cols else []
        )

    def load_normalized_data(
        self,
        adata: AnnData,
        compound_id_obs_col: str,
        layer: Optional[str] = None,
        sample_id_obs_cols: Optional[List[str]] = None,
        signature_id_obs_cols: Optional[List[str]] = None,
    ):
        if layer and layer in adata.layers:
            temp_adata = adata.copy()
            temp_adata.X = temp_adata.layers[layer]
        else:
            temp_adata = adata.copy()

        if sample_id_obs_cols is None:
            sample_id_obs_cols = []
        if signature_id_obs_cols is None:
            signature_id_obs_cols = []

        sample_pseudobulk_cols = (
            [compound_id_obs_col] + sample_id_obs_cols + signature_id_obs_cols
        )
        pb_adata = pseudobulk_adata(
            temp_adata,
            sample_id_obs_cols=sample_pseudobulk_cols,
            method="mean",
        )

        self.expr = pb_adata
        self._compound_id_obs_col = compound_id_obs_col
        self._sample_id_obs_cols = sample_id_obs_cols if sample_id_obs_cols else []
        self._signature_id_obs_cols = (
            signature_id_obs_cols if signature_id_obs_cols else []
        )

    def load_phenotypic_readouts(
        self,
        readouts: Union[pd.DataFrame, pd.Series],
        readout_col: Optional[str] = None,
        compound_id_col: Optional[str] = None,
    ):
        if isinstance(readouts, pd.Series):
            readout_series = readouts.copy()
        elif isinstance(readouts, pd.DataFrame):
            if readout_col is None:
                raise ValueError(
                    "readout_col must be specified when readouts is a DataFrame"
                )

            if compound_id_col is None:
                compound_ids = readouts.index
            else:
                compound_ids = readouts[compound_id_col]

            readout_series = pd.Series(readouts[readout_col].values, index=compound_ids)
        else:
            raise ValueError("readouts must be a pandas DataFrame or Series")

        self.readouts = readout_series.groupby(readout_series.index).mean()
        return self.readouts

    def compute_learned_signatures(self, corr_method: str = "pearson"):
        if self.expr is None:
            raise ValueError("Expression data not loaded. Call load_counts_data first.")
        if self.readouts is None:
            raise ValueError(
                "Phenotypic readouts not loaded. Call load_phenotypic_readouts first."
            )
        if corr_method != "pearson":
            raise ValueError("Only pearson correlation is currently supported.")

        group_cols = [self._compound_id_obs_col] + self._signature_id_obs_cols
        expr_obs = self.expr.obs.copy()
        expr_obs["_group_key"] = expr_obs[group_cols].astype(str).agg("|".join, axis=1)

        learned_rows = []
        learned_obs = []

        for signature_key, subset_idx in expr_obs.groupby("_group_key").groups.items():
            subset = self.expr[list(subset_idx)].copy()
            compound_ids = subset.obs[self._compound_id_obs_col].astype(str)
            matching = compound_ids.isin(self.readouts.index.astype(str))
            subset = subset[matching].copy()
            compound_ids = subset.obs[self._compound_id_obs_col].astype(str)

            if subset.n_obs < 2:
                continue

            X = subset.X.toarray() if sp.issparse(subset.X) else np.asarray(subset.X)
            y = self.readouts.reindex(compound_ids).to_numpy()

            corrs = np.zeros(subset.n_vars, dtype=float)
            for gene_idx in range(subset.n_vars):
                corrs[gene_idx] = pearsonr(X[:, gene_idx], y).statistic

            learned_rows.append(corrs)
            learned_obs.append(signature_key)

        if not learned_rows:
            raise ValueError("Unable to compute learned signatures from the provided data.")

        self.learned_signatures = AnnData(
            X=np.vstack(learned_rows),
            obs=pd.DataFrame(index=learned_obs),
            var=pd.DataFrame(index=self.expr.var_names),
        )
        return self.learned_signatures

    def compute_refined_signatures(
        self,
        learning_rate: float = 0.5,
        scale_learned_sig: bool = True,
    ):
        if self.learned_signatures is None:
            raise ValueError(
                "Learned signatures not computed. Call compute_learned_signatures first."
            )

        learned = self.learned_signatures.X.copy()
        starting = self.starting_signature.reindex(
            self.learned_signatures.var_names, fill_value=0.0
        ).to_numpy()

        if scale_learned_sig:
            starting_std = np.std(starting)
            learned_std = np.std(learned, axis=1, keepdims=True)
            learned_std[learned_std == 0] = 1.0
            learned = learned * (starting_std / learned_std)

        refined = (1 - learning_rate) * starting[np.newaxis, :] + learning_rate * learned

        self.refined_signatures = AnnData(
            X=refined,
            obs=self.learned_signatures.obs.copy(),
            var=self.learned_signatures.var.copy(),
        )
        return self.refined_signatures
