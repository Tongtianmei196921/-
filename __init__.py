"""Compatibility exports for the source-root package layout."""

from .drugreflector import (
    DrugReflector,
    EnsembleModel,
    SignatureRefinement,
    compute_vscore_two_groups,
    compute_vscores,
    compute_vscores_adata,
    create_synthetic_gene_expression,
    load_h5ad_file,
    nnFC,
    pseudobulk_adata,
)

__version__ = "1.0.0"
__all__ = [
    "DrugReflector",
    "EnsembleModel",
    "nnFC",
    "load_h5ad_file",
    "create_synthetic_gene_expression",
    "compute_vscores",
    "compute_vscore_two_groups",
    "compute_vscores_adata",
    "pseudobulk_adata",
    "SignatureRefinement",
]
