"""Helpers for fetching, parsing, and analyzing GEO series matrices."""

from __future__ import annotations

import csv
import gzip
import io
import os
import re
from http.client import IncompleteRead
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
from urllib.request import urlopen

import numpy as np
import pandas as pd
from anndata import AnnData

from .utils import compute_vscores_adata

APP_DIR = Path(__file__).resolve().parent
REPO_ROOT = APP_DIR.parent
GEO_CACHE_DIR = REPO_ROOT / ".geo_cache"
MAX_GEO_DOWNLOAD_BYTES = int(os.getenv("DRUGREFLECTOR_MAX_GEO_DOWNLOAD_MB", "500")) * 1024 * 1024

_CONTROL_HINTS = (
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
_CASE_HINTS = (
    "tumor",
    "cancer",
    "disease",
    "treated",
    "resistant",
    "mutant",
    "knockout",
    "ko",
    "kd",
    "disease state",
    "patient",
    "infected",
    "stimulated",
    "drug",
)
_GROUP_COLUMN_HINTS = (
    "group",
    "condition",
    "phenotype",
    "source_name",
    "source",
    "characteristics",
    "title",
    "disease",
    "status",
    "treatment",
    "sample_type",
)
_SYMBOL_COLUMN_HINTS = (
    "gene symbol",
    "gene_symbol",
    "symbol",
    "gene symbol ch1",
    "genesymbol",
    "gene",
)
_INVALID_SYMBOLS = {"", "---", "na", "n/a", "null", "nan", "none"}


@dataclass
class GeoDataset:
    accession: str
    sample_metadata: pd.DataFrame
    expression_by_probe: pd.DataFrame
    expression_by_gene: pd.DataFrame
    platform_id: str | None
    symbol_source: str
    used_log2: bool
    organism: str | None
    expression_source: str
    ortholog_mapping: dict[str, object] | None


@dataclass
class GeoGrouping:
    group_column: str
    group1_value: str
    group2_value: str
    group1_count: int
    group2_count: int
    mode: str


def _bucket_prefix(accession: str, prefix: str) -> str:
    digits = accession[len(prefix) :]
    if not digits.isdigit():
        raise ValueError(f"Unsupported GEO accession: {accession}")
    head = digits[:-3]
    return f"{prefix}{head}nnn" if head else f"{prefix}nnn"


def _series_matrix_url(accession: str) -> str:
    accession = accession.upper()
    return (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{_bucket_prefix(accession, 'GSE')}/{accession}/matrix/{accession}_series_matrix.txt.gz"
    )


def _platform_annotation_url(platform_id: str) -> str:
    platform_id = platform_id.upper()
    return (
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/"
        f"{_bucket_prefix(platform_id, 'GPL')}/{platform_id}/annot/{platform_id}.annot.gz"
    )


def _series_supplementary_url(accession: str) -> str:
    accession = accession.upper()
    return (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{_bucket_prefix(accession, 'GSE')}/{accession}/suppl/"
    )


def _normalize_geo_url(url: str) -> str:
    cleaned = str(url).strip().strip('"')
    if cleaned.lower() in _INVALID_SYMBOLS:
        return ""
    if not re.match(r"^(https?|ftp)://", cleaned, flags=re.IGNORECASE):
        return ""
    if cleaned.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
        return "https://ftp.ncbi.nlm.nih.gov/" + cleaned.removeprefix("ftp://ftp.ncbi.nlm.nih.gov/")
    return cleaned


def _cached_download(url: str, cache_name: str) -> bytes:
    GEO_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = GEO_CACHE_DIR / cache_name
    if cache_path.exists():
        return cache_path.read_bytes()

    last_error: Exception | None = None
    for _ in range(3):
        try:
            with urlopen(url, timeout=120) as response:
                data = response.read(MAX_GEO_DOWNLOAD_BYTES + 1)
            if len(data) > MAX_GEO_DOWNLOAD_BYTES:
                limit_mb = MAX_GEO_DOWNLOAD_BYTES // (1024 * 1024)
                raise ValueError(f"GEO resource is too large to process safely: {url} (>{limit_mb} MB)")
            cache_path.write_bytes(data)
            return data
        except IncompleteRead as exc:  # pragma: no cover - network dependent
            last_error = exc
            continue
        except HTTPError as exc:  # pragma: no cover - network dependent
            raise ValueError(f"Could not download GEO resource: {url} ({exc.code})") from exc
        except URLError as exc:  # pragma: no cover - network dependent
            raise ValueError(f"Could not reach GEO resource: {url}") from exc

    if last_error is not None:  # pragma: no cover - network dependent
        raise ValueError(f"Could not fully download GEO resource after retries: {url}") from last_error

    raise ValueError(f"Could not download GEO resource: {url}")


def _gunzip_text(raw: bytes) -> str:
    return gzip.decompress(raw).decode("utf-8", errors="replace")


def _parse_tsv_line(line: str) -> list[str]:
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def _parse_series_matrix(accession: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = _cached_download(
        _series_matrix_url(accession),
        f"{accession.upper()}_series_matrix.txt.gz",
    )
    text = _gunzip_text(raw)

    metadata_rows: dict[str, list[str]] = {}
    key_counts: dict[str, int] = {}
    table_lines: list[str] = []
    in_table = False

    for line in text.splitlines():
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if in_table:
            table_lines.append(line)
            continue
        if line.startswith("!Sample_"):
            values = _parse_tsv_line(line)
            key = values[0][len("!Sample_") :]
            key_counts[key] = key_counts.get(key, 0) + 1
            if key_counts[key] > 1:
                key = f"{key}_{key_counts[key]}"
            metadata_rows[key] = values[1:]

    if not table_lines:
        raise ValueError(f"GEO accession {accession.upper()} does not contain a series matrix table.")

    matrix = pd.read_csv(io.StringIO("\n".join(table_lines)), sep="\t")
    matrix = matrix.rename(columns={matrix.columns[0]: "ID_REF"})
    matrix["ID_REF"] = matrix["ID_REF"].astype(str).str.strip().str.replace('"', "")
    matrix = matrix.set_index("ID_REF")
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.dropna(axis=0, how="all")

    if not metadata_rows:
        raise ValueError(f"GEO accession {accession.upper()} does not expose sample metadata.")

    sample_ids = metadata_rows.get("geo_accession", list(matrix.columns))
    obs = pd.DataFrame(metadata_rows, index=sample_ids)
    obs.index.name = "sample_id"

    if len(obs.index) != len(matrix.columns):
        obs = obs.reindex(matrix.columns)
        obs.index.name = "sample_id"

    matrix = matrix.loc[:, obs.index.tolist()]
    return obs, matrix


def _parse_platform_annotation(platform_id: str) -> pd.DataFrame:
    raw = _cached_download(
        _platform_annotation_url(platform_id),
        f"{platform_id.upper()}.annot.gz",
    )
    text = _gunzip_text(raw)

    table_lines: list[str] = []
    in_table = False
    for line in text.splitlines():
        if line.startswith("!platform_table_begin"):
            in_table = True
            continue
        if line.startswith("!platform_table_end"):
            break
        if in_table:
            table_lines.append(line)

    if not table_lines:
        raise ValueError(f"Platform annotation for {platform_id.upper()} does not contain a table.")

    annot = pd.read_csv(io.StringIO("\n".join(table_lines)), sep="\t")
    annot = annot.rename(columns={annot.columns[0]: "ID"})
    annot["ID"] = annot["ID"].astype(str).str.strip().str.replace('"', "")
    return annot


def _canonical_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.strip().lower()).strip()


def _extract_gene_symbol(raw_value: object) -> str | None:
    if raw_value is None or pd.isna(raw_value):
        return None

    text = str(raw_value).strip().replace('"', "")
    if not text:
        return None

    for token in re.split(r"///|//|;|,|\|", text):
        cleaned = token.strip().upper()
        cleaned = re.sub(r"\s+", "", cleaned)
        if not cleaned:
            continue
        if cleaned.lower() in _INVALID_SYMBOLS:
            continue
        if not re.fullmatch(r"[A-Z0-9\-]{2,20}", cleaned):
            continue
        if cleaned.startswith(("AFFX", "ILMN_", "A_", "ENSG", "ENST", "XM_", "NM_")):
            continue
        if cleaned.isdigit():
            continue
        return cleaned
    return None


def _looks_like_gene_symbols(index: pd.Index) -> bool:
    values = pd.Index(index.astype(str).str.strip()).unique()
    if len(values) == 0:
        return False

    values = values[: min(len(values), 500)]
    valid = 0
    for value in values:
        upper = value.upper()
        if upper.startswith(("AFFX", "ILMN_", "A_", "NM_", "XM_")):
            continue
        if re.fullmatch(r"ENS[A-Z]*[GT]\d+(?:\.\d+)?", upper):
            continue
        if "_" in upper:
            continue
        if re.fullmatch(r"[A-Z0-9\-]{2,20}", upper) and not upper.isdigit():
            valid += 1
    return valid / max(len(values), 1) >= 0.6


def _looks_like_ensembl_ids(index: pd.Index) -> bool:
    values = pd.Index(index.astype(str).str.strip()).unique()
    if len(values) == 0:
        return False
    values = values[: min(len(values), 500)]
    matches = 0
    for value in values:
        if re.fullmatch(r"ENS[A-Z]*G\d+(?:\.\d+)?", value.upper()):
            matches += 1
    return matches / max(len(values), 1) >= 0.6


def _select_symbol_column(columns: pd.Index) -> str | None:
    normalized = {_canonical_name(str(column)): str(column) for column in columns}
    for hint in _SYMBOL_COLUMN_HINTS:
        if hint in normalized:
            return normalized[hint]
    for column in columns:
        name = _canonical_name(str(column))
        if "symbol" in name and "gene" in name:
            return str(column)
    return None


def _collapse_gene_expression(frame: pd.DataFrame, symbols: pd.Series) -> pd.DataFrame:
    valid = symbols.dropna()
    valid = valid[valid.index.isin(frame.index)]
    if valid.empty:
        raise ValueError("No gene symbols could be recovered from the GEO platform annotation.")

    subset = frame.loc[valid.index].copy()
    subset.index = valid.loc[subset.index]
    subset = subset.groupby(level=0).mean()
    subset = subset[~subset.index.duplicated(keep="first")]
    return subset.sort_index()


def _log2_if_needed(frame: pd.DataFrame) -> tuple[pd.DataFrame, bool]:
    values = frame.to_numpy(dtype=float, copy=False)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise ValueError("GEO expression matrix contains no finite values.")
    if finite.min() < 0:
        return frame, False
    if np.nanpercentile(finite, 99) > 100 or np.nanmax(finite) > 1000:
        return np.log2(frame + 1.0), True
    return frame, False


def _detect_organism(sample_metadata: pd.DataFrame) -> str | None:
    organism_columns = [column for column in sample_metadata.columns if "organism" in _canonical_name(column)]
    for column in organism_columns:
        values = [
            value
            for value in sample_metadata[column].fillna("").astype(str).str.strip().unique().tolist()
            if value
        ]
        if len(values) == 1:
            return values[0]
    return None


def _pick_supplementary_column(sample_metadata: pd.DataFrame) -> str | None:
    supplementary_columns = [
        column
        for column in sample_metadata.columns
        if "supplementary file" in _canonical_name(column)
    ]
    ranked: list[tuple[tuple[int, int, int], str]] = []
    for column in supplementary_columns:
        values = sample_metadata[column].fillna("").astype(str).str.strip()
        non_empty = [_normalize_geo_url(value) for value in values.tolist()]
        non_empty = [value for value in non_empty if value]
        if len(non_empty) != len(sample_metadata.index):
            continue
        joined = " ".join(non_empty).lower()
        rank = (
            1 if "normalized" in joined else 0,
            1 if "count" in joined else 0,
            len(non_empty),
        )
        ranked.append((rank, column))
    ranked.sort(reverse=True)
    return ranked[0][1] if ranked else None


def _parse_supplementary_expression_file(raw: bytes, sample_id: str) -> pd.Series:
    if raw[:2] == b"\x1f\x8b":
        text = gzip.decompress(raw).decode("utf-8", errors="replace")
    else:
        text = raw.decode("utf-8", errors="replace")

    frame = pd.read_csv(
        io.StringIO(text),
        sep=r"\s+|,",
        engine="python",
        comment="#",
        header=None,
    )
    if frame.shape[1] < 2:
        raise ValueError(f"Supplementary expression file for {sample_id} does not contain at least two columns.")

    frame = frame.iloc[:, :2].copy()
    frame.columns = ["gene_id", "value"]
    frame["gene_id"] = frame["gene_id"].astype(str).str.strip().str.replace('"', "")
    frame["value"] = pd.to_numeric(frame["value"], errors="coerce")
    frame = frame.dropna(subset=["gene_id", "value"])
    frame = frame[frame["gene_id"] != ""]
    if frame.empty:
        raise ValueError(f"Supplementary expression file for {sample_id} is empty after parsing.")
    return pd.Series(frame["value"].to_numpy(dtype=float), index=frame["gene_id"], name=sample_id)


def _load_supplementary_expression(sample_metadata: pd.DataFrame) -> tuple[pd.DataFrame, str] | None:
    column = _pick_supplementary_column(sample_metadata)
    if not column:
        return None

    series_by_sample: dict[str, pd.Series] = {}
    for sample_id, url in sample_metadata[column].fillna("").astype(str).items():
        normalized_url = _normalize_geo_url(url)
        if not normalized_url:
            return None
        cache_name = Path(normalized_url).name or f"{sample_id}.txt"
        raw = _cached_download(normalized_url, cache_name)
        series_by_sample[str(sample_id)] = _parse_supplementary_expression_file(raw, str(sample_id))

    frame = pd.DataFrame(series_by_sample)
    frame = frame.dropna(axis=0, how="all").dropna(axis=1, how="all")
    if frame.empty:
        raise ValueError("Could not build an expression matrix from GEO supplementary files.")
    return frame, column


def _list_series_supplementary_files(accession: str) -> list[str]:
    base_url = _series_supplementary_url(accession)
    raw = _cached_download(base_url, f"{accession.upper()}_suppl_index.html")
    text = raw.decode("utf-8", errors="replace")
    files: list[str] = []
    for href in re.findall(r'href="([^"]+)"', text, flags=re.IGNORECASE):
        if href.startswith("?") or href.startswith("/") or href == "../":
            continue
        normalized = _normalize_geo_url(urljoin(base_url, href))
        if normalized:
            files.append(normalized)
    return files


def _read_supplementary_table(url: str, accession: str) -> pd.DataFrame:
    cache_name = f"{accession.upper()}_{Path(url).name}"
    raw = _cached_download(url, cache_name)
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    lower = Path(url).name.lower()
    if lower.endswith((".csv", ".csv.gz")):
        return pd.read_csv(io.BytesIO(raw))
    if lower.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        return pd.read_csv(io.BytesIO(raw), sep=None, engine="python")
    raise ValueError(f"Unsupported supplementary table format: {Path(url).name}")


def _select_differential_table(accession: str) -> tuple[str, pd.DataFrame] | None:
    candidates = []
    for url in _list_series_supplementary_files(accession):
        name = Path(url).name.lower()
        if not name.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz")):
            continue
        score = 0
        if any(token in name for token in ["regression", "differential", "diff", "limma", "deg", "result"]):
            score += 10
        if any(token in name for token in ["count", "cnt", "norm"]):
            score -= 4
        candidates.append((score, url))
    for _, url in sorted(candidates, reverse=True):
        try:
            frame = _read_supplementary_table(url, accession)
        except Exception:
            continue
        columns = [_canonical_name(column) for column in frame.columns]
        has_feature = any(name in columns for name in ["gene", "genesymbol", "symbol", "mirna", "feature"])
        has_score = any("logfc" in name or name in {"stat", "score", "t"} for name in columns)
        if has_feature and has_score:
            return url, frame
    return None


def _feature_column(frame: pd.DataFrame) -> str | None:
    preferred = ["gene", "genesymbol", "symbol", "geneid", "mirna", "feature"]
    normalized = {_canonical_name(column): str(column) for column in frame.columns}
    for name in preferred:
        if name in normalized:
            return normalized[name]
    return str(frame.columns[0]) if len(frame.columns) else None


def _score_column(frame: pd.DataFrame) -> str | None:
    columns = [str(column) for column in frame.columns]
    preferred_tokens = ["ageadj.logfc", "logfc", "stat", "score", "t"]
    ranked: list[tuple[int, str]] = []
    for column in columns:
        name = _canonical_name(column)
        if "pvalue" in name or "qvalue" in name or "padj" in name:
            continue
        for idx, token in enumerate(preferred_tokens):
            if token in name:
                ranked.append((len(preferred_tokens) - idx, column))
                break
    ranked.sort(reverse=True)
    return ranked[0][1] if ranked else None


def _looks_like_mirna_features(values: pd.Series, column: str) -> bool:
    column_name = _canonical_name(column)
    if "mirna" in column_name:
        return True
    sample = values.dropna().astype(str).str.lower().head(50)
    if sample.empty:
        return False
    return float(sample.str.match(r"^(hsa-)?(mir|miR|let-|hsa-let-)").mean()) >= 0.5


def _load_series_differential_signature(accession: str) -> tuple[pd.DataFrame, dict[str, object]]:
    selected = _select_differential_table(accession)
    if selected is None:
        raise ValueError(
            f"GEO accession {accession} does not expose an expression matrix in the series matrix "
            "and no per-sample supplementary expression files or gene-level differential table were found."
        )
    url, table = selected
    feature_column = _feature_column(table)
    score_column = _score_column(table)
    if not feature_column or not score_column:
        raise ValueError(f"Could not identify feature and score columns in {Path(url).name}.")
    if _looks_like_mirna_features(table[feature_column], feature_column):
        raise ValueError(
            f"GEO accession {accession} provides a differential miRNA table ({Path(url).name}), "
            "not a human gene-level differential expression table. DrugReflector requires gene symbols; "
            "miRNA IDs cannot be objectively treated as genes without an additional validated target-mapping step."
        )

    work = table[[feature_column, score_column]].copy()
    work.columns = ["gene", "score"]
    work["gene"] = work["gene"].map(_extract_gene_symbol)
    work["score"] = pd.to_numeric(work["score"], errors="coerce")
    work = work.dropna(subset=["gene", "score"])
    work = work.groupby("gene", sort=False)["score"].mean()
    if work.empty:
        raise ValueError(f"Differential table {Path(url).name} did not contain usable gene-level scores.")
    frame = pd.DataFrame([work], index=[f"{accession.upper()}:{Path(url).name}:{score_column}"])
    metadata = {
        "accession": accession.upper(),
        "platform_id": None,
        "organism": "Homo sapiens",
        "symbol_source": f"supplementary_differential:{feature_column}",
        "expression_source": f"supplementary_differential:{Path(url).name}",
        "ortholog_mapping": None,
        "used_log2": False,
        "group_column": None,
        "group1_value": None,
        "group2_value": None,
        "group1_count": None,
        "group2_count": None,
        "n_samples_total": None,
        "n_genes": int(frame.shape[1]),
        "mode": "author_differential_table",
        "differential_score_column": score_column,
        "differential_table_url": url,
    }
    return frame, metadata


def _map_mouse_ensembl_to_human_symbols(
    gene_ids: pd.Index,
    *,
    cache_key: str,
) -> tuple[pd.Series, dict[str, object]]:
    gene_list = pd.Index(gene_ids.astype(str).str.strip().str.upper()).unique()
    gene_list = gene_list[gene_list.str.fullmatch(r"ENSMUSG\d+(?:\.\d+)?", na=False)]
    if len(gene_list) == 0:
        return pd.Series(dtype=object), {
            "source": "ensembl_biomart_mouse_to_human",
            "input_genes": 0,
            "mapped_input_genes": 0,
            "unique_human_symbols": 0,
            "orthology_types": [],
        }

    ensembl_raw = _cached_download(
        "https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt",
        "MGI_MRK_ENSEMBL.rpt",
    )
    homology_raw = _cached_download(
        "https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt",
        "MGI_HOM_ProteinCoding.rpt",
    )
    ensembl = pd.read_csv(
        io.StringIO(ensembl_raw.decode("utf-8", errors="replace")),
        sep="\t",
        header=None,
        usecols=[0, 5],
        names=["mgi_id", "mouse_ensembl_gene_id"],
    )
    homology = pd.read_csv(
        io.StringIO(homology_raw.decode("utf-8", errors="replace")),
        sep="\t",
        header=None,
        usecols=[0, 4],
        names=["mgi_id", "human_gene_name"],
    )
    mapping = ensembl.merge(homology, on="mgi_id", how="inner")
    mapping["mouse_ensembl_gene_id"] = mapping["mouse_ensembl_gene_id"].astype(str).str.strip().str.upper()
    mapping["human_gene_name"] = mapping["human_gene_name"].map(_extract_gene_symbol)
    mapping = mapping[mapping["human_gene_name"].notna()].copy()
    mapping = mapping[mapping["mouse_ensembl_gene_id"].isin(gene_list)].copy()
    mapping = mapping.drop_duplicates(subset=["mouse_ensembl_gene_id"], keep="first")

    if mapping.empty:
        return pd.Series(dtype=object), {
            "source": "mgi_mouse_human_one_to_one",
            "input_genes": int(len(gene_list)),
            "mapped_input_genes": 0,
            "unique_human_symbols": 0,
            "orthology_types": ["one_to_one"],
        }

    symbol_map = pd.Series(
        mapping["human_gene_name"].to_numpy(),
        index=mapping["mouse_ensembl_gene_id"],
    )
    stats = {
        "source": "mgi_mouse_human_one_to_one",
        "input_genes": int(len(gene_list)),
        "mapped_input_genes": int(mapping["mouse_ensembl_gene_id"].nunique()),
        "unique_human_symbols": int(mapping["human_gene_name"].nunique()),
        "orthology_types": ["one_to_one"],
    }
    return symbol_map, stats


def load_geo_dataset(accession: str) -> GeoDataset:
    accession = accession.strip().upper()
    if not re.fullmatch(r"GSE\d+", accession):
        raise ValueError("Please provide a GEO series accession such as GSE6631.")

    sample_metadata, expression_by_probe = _parse_series_matrix(accession)
    platform_ids = sample_metadata.get("platform_id", pd.Series(dtype=object)).dropna().unique()
    platform_id = str(platform_ids[0]) if len(platform_ids) == 1 else None
    organism = _detect_organism(sample_metadata)
    ortholog_mapping: dict[str, object] | None = None

    expression_source = "series_matrix"
    if expression_by_probe.empty:
        supplementary = _load_supplementary_expression(sample_metadata)
        if supplementary is None:
            raise ValueError(
                f"GEO accession {accession} does not expose an expression matrix in the series matrix "
                "and no per-sample supplementary expression files were found."
            )
        expression_by_probe, supplementary_column = supplementary
        expression_source = f"supplementary:{supplementary_column}"

    if _looks_like_gene_symbols(expression_by_probe.index):
        symbols = pd.Series(
            expression_by_probe.index.astype(str).str.upper(),
            index=expression_by_probe.index,
        )
        expression_by_gene = _collapse_gene_expression(expression_by_probe, symbols)
        symbol_source = "series_matrix"
    elif _looks_like_ensembl_ids(expression_by_probe.index):
        expression_by_probe.index = pd.Index(expression_by_probe.index.astype(str).str.upper())
        if organism and organism.lower() == "mus musculus":
            symbols, ortholog_mapping = _map_mouse_ensembl_to_human_symbols(
                expression_by_probe.index,
                cache_key=accession,
            )
            expression_by_gene = _collapse_gene_expression(expression_by_probe, symbols)
            symbol_source = "mouse_to_human_orthologs"
        else:
            expression_by_gene = expression_by_probe.copy()
            symbol_source = "ensembl_ids"
    else:
        if not platform_id:
            raise ValueError(
                "This GEO dataset appears to use probe IDs and does not expose a single platform ID."
            )
        annotation = _parse_platform_annotation(platform_id)
        symbol_column = _select_symbol_column(annotation.columns)
        if not symbol_column:
            raise ValueError(f"Could not find a gene symbol column in platform {platform_id}.")
        symbols = (
            annotation.set_index("ID")[symbol_column]
            .map(_extract_gene_symbol)
            .dropna()
        )
        expression_by_gene = _collapse_gene_expression(expression_by_probe, symbols)
        symbol_source = f"{platform_id}:{symbol_column}"

    expression_by_gene, used_log2 = _log2_if_needed(expression_by_gene)
    return GeoDataset(
        accession=accession,
        sample_metadata=sample_metadata,
        expression_by_probe=expression_by_probe,
        expression_by_gene=expression_by_gene,
        platform_id=platform_id,
        symbol_source=symbol_source,
        used_log2=used_log2,
        organism=organism,
        expression_source=expression_source,
        ortholog_mapping=ortholog_mapping,
    )


def _score_group_column(column: str, values: pd.Series) -> float:
    non_null = values.fillna("").astype(str).str.strip()
    unique = [value for value in non_null.unique().tolist() if value]
    unique_count = len(unique)
    if unique_count < 2 or unique_count > min(8, max(len(values) - 1, 2)):
        return -1.0

    score = 0.0
    if unique_count == 2:
        score += 40.0
    elif unique_count <= 4:
        score += 10.0

    normalized_name = _canonical_name(column)
    for hint in _GROUP_COLUMN_HINTS:
        if hint in normalized_name:
            score += 8.0

    joined_values = " ".join(value.lower() for value in unique)
    if any(hint in joined_values for hint in _CONTROL_HINTS):
        score += 10.0
    if any(hint in joined_values for hint in _CASE_HINTS):
        score += 10.0

    counts = non_null.value_counts()
    if counts.min() >= 2:
        score += 4.0
    return score


def _order_groups(values: list[str]) -> tuple[str, str]:
    def score(value: str) -> tuple[int, int, str]:
        lower = value.lower()
        control = sum(hint in lower for hint in _CONTROL_HINTS)
        case = sum(hint in lower for hint in _CASE_HINTS)
        # Lower sort key should become group1/reference.
        return (0 if control >= case else 1, -control, lower)

    ordered = sorted(values, key=score)
    return ordered[0], ordered[1]


def infer_geo_grouping(
    sample_metadata: pd.DataFrame,
    *,
    group_column: str | None = None,
    group1_value: str | None = None,
    group2_value: str | None = None,
    control_keyword: str | None = None,
    case_keyword: str | None = None,
) -> GeoGrouping:
    obs = sample_metadata.copy()

    if group_column and group1_value and group2_value:
        if group_column not in obs.columns:
            raise ValueError(f"GEO metadata column '{group_column}' was not found.")
        group_values = obs[group_column].fillna("").astype(str)
        count1 = int((group_values == group1_value).sum())
        count2 = int((group_values == group2_value).sum())
        if count1 == 0 or count2 == 0:
            raise ValueError("Selected GEO group values do not match any samples.")
        return GeoGrouping(group_column, group1_value, group2_value, count1, count2, "manual")

    if control_keyword and case_keyword:
        haystack = obs.fillna("").astype(str).agg(" | ".join, axis=1).str.lower()
        control_keyword = control_keyword.strip().lower()
        case_keyword = case_keyword.strip().lower()
        control_mask = haystack.str.contains(re.escape(control_keyword), regex=True)
        case_mask = haystack.str.contains(re.escape(case_keyword), regex=True)
        overlap = control_mask & case_mask
        control_mask = control_mask & ~overlap
        case_mask = case_mask & ~overlap
        if not control_mask.any() or not case_mask.any():
            raise ValueError(
                "Could not match the provided GEO keywords to two non-empty sample groups."
            )
        return GeoGrouping(
            "_keyword_match",
            control_keyword,
            case_keyword,
            int(control_mask.sum()),
            int(case_mask.sum()),
            "keyword",
        )

    ranked: list[tuple[float, str]] = []
    for column in obs.columns:
        ranked.append((_score_group_column(column, obs[column]), column))
    ranked.sort(reverse=True)
    if not ranked or ranked[0][0] < 0:
        raise ValueError(
            "Could not auto-detect a GEO grouping column. Please provide group keywords."
        )

    best_column = ranked[0][1]
    values = [
        value
        for value in obs[best_column].fillna("").astype(str).str.strip().unique().tolist()
        if value
    ]
    if len(values) != 2:
        raise ValueError(
            f"Auto-detected GEO column '{best_column}' has {len(values)} groups; "
            "please specify group keywords or explicit values."
        )

    group1_value, group2_value = _order_groups(values)
    group_values = obs[best_column].fillna("").astype(str)
    return GeoGrouping(
        best_column,
        group1_value,
        group2_value,
        int((group_values == group1_value).sum()),
        int((group_values == group2_value).sum()),
        "auto",
    )


def preview_geo(accession: str) -> dict[str, object]:
    accession = accession.strip().upper()
    try:
        dataset = load_geo_dataset(accession)
    except ValueError as exc:
        if "does not expose an expression matrix" not in str(exc):
            raise
        frame, metadata = _load_series_differential_signature(accession)
        return {
            "accession": accession,
            "n_samples": 0,
            "n_probe_rows": int(frame.shape[1]),
            "n_genes": int(frame.shape[1]),
            "platform_id": metadata.get("platform_id"),
            "organism": metadata.get("organism"),
            "symbol_source": metadata.get("symbol_source"),
            "expression_source": metadata.get("expression_source"),
            "ortholog_mapping": metadata.get("ortholog_mapping"),
            "used_log2": metadata.get("used_log2"),
            "detected_grouping": None,
            "candidate_columns": [],
            "sample_preview": [],
            "mode": metadata.get("mode"),
            "differential_score_column": metadata.get("differential_score_column"),
            "differential_table_url": metadata.get("differential_table_url"),
        }

    candidate_columns = []
    for column in dataset.sample_metadata.columns:
        score = _score_group_column(column, dataset.sample_metadata[column])
        if score < 0:
            continue
        values = [
            value
            for value in dataset.sample_metadata[column].fillna("").astype(str).str.strip().unique().tolist()
            if value
        ]
        candidate_columns.append(
            {"column": column, "n_unique": len(values), "values": values[:6], "score": round(score, 2)}
        )

    candidate_columns.sort(key=lambda item: item["score"], reverse=True)

    detected = None
    try:
        grouping = infer_geo_grouping(dataset.sample_metadata)
        detected = {
            "group_column": grouping.group_column,
            "group1_value": grouping.group1_value,
            "group2_value": grouping.group2_value,
            "group1_count": grouping.group1_count,
            "group2_count": grouping.group2_count,
            "mode": grouping.mode,
        }
    except ValueError:
        detected = None

    preview_rows = dataset.sample_metadata.head(8).reset_index().fillna("").to_dict(orient="records")
    return {
        "accession": dataset.accession,
        "n_samples": int(dataset.sample_metadata.shape[0]),
        "n_probe_rows": int(dataset.expression_by_probe.shape[0]),
        "n_genes": int(dataset.expression_by_gene.shape[0]),
        "platform_id": dataset.platform_id,
        "organism": dataset.organism,
        "symbol_source": dataset.symbol_source,
        "expression_source": dataset.expression_source,
        "ortholog_mapping": dataset.ortholog_mapping,
        "used_log2": dataset.used_log2,
        "detected_grouping": detected,
        "candidate_columns": candidate_columns[:8],
        "sample_preview": preview_rows,
    }


def build_geo_signature(
    accession: str,
    *,
    group_column: str | None = None,
    group1_value: str | None = None,
    group2_value: str | None = None,
    control_keyword: str | None = None,
    case_keyword: str | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    accession = accession.strip().upper()
    try:
        dataset = load_geo_dataset(accession)
    except ValueError as exc:
        if "does not expose an expression matrix" in str(exc):
            return _load_series_differential_signature(accession)
        raise
    if dataset.organism and dataset.organism.lower() != "homo sapiens" and dataset.symbol_source != "mouse_to_human_orthologs":
        raise ValueError(
            f"GEO accession {dataset.accession} is {dataset.organism}. "
            "The current automatic DrugReflector GEO pipeline only supports human datasets "
            "unless genes have already been mapped to human symbols objectively."
        )
    if dataset.symbol_source == "ensembl_ids":
        raise ValueError(
            f"GEO accession {dataset.accession} exposes Ensembl-style gene IDs. "
            "Automatic GEO import currently requires human gene symbols after parsing; "
            "please pre-map genes to HGNC symbols or upload a prepared signature file."
        )
    grouping = infer_geo_grouping(
        dataset.sample_metadata,
        group_column=group_column,
        group1_value=group1_value,
        group2_value=group2_value,
        control_keyword=control_keyword,
        case_keyword=case_keyword,
    )

    obs = dataset.sample_metadata.copy()
    if grouping.mode == "keyword":
        haystack = obs.fillna("").astype(str).agg(" | ".join, axis=1).str.lower()
        obs["_keyword_match"] = np.where(
            haystack.str.contains(re.escape(grouping.group1_value), regex=True),
            grouping.group1_value,
            np.where(
                haystack.str.contains(re.escape(grouping.group2_value), regex=True),
                grouping.group2_value,
                "",
            ),
        )
        obs = obs.loc[obs["_keyword_match"] != ""].copy()
        expression = dataset.expression_by_gene.loc[:, obs.index]
    else:
        expression = dataset.expression_by_gene

    adata = AnnData(
        X=expression.T.to_numpy(dtype=np.float32, copy=True),
        obs=obs.copy(),
        var=pd.DataFrame(index=expression.index.astype(str)),
    )
    adata.obs.index = obs.index.astype(str)

    vscores = compute_vscores_adata(
        adata,
        grouping.group_column,
        grouping.group1_value,
        grouping.group2_value,
    )
    signature_name = f"{dataset.accession}:{grouping.group1_value}->{grouping.group2_value}"
    frame = pd.DataFrame([vscores], index=[signature_name])

    metadata = {
        "accession": dataset.accession,
        "platform_id": dataset.platform_id,
        "organism": dataset.organism,
        "symbol_source": dataset.symbol_source,
        "expression_source": dataset.expression_source,
        "ortholog_mapping": dataset.ortholog_mapping,
        "used_log2": dataset.used_log2,
        "group_column": grouping.group_column,
        "group1_value": grouping.group1_value,
        "group2_value": grouping.group2_value,
        "group1_count": grouping.group1_count,
        "group2_count": grouping.group2_count,
        "n_samples_total": int(dataset.sample_metadata.shape[0]),
        "n_genes": int(frame.shape[1]),
        "mode": grouping.mode,
    }
    return frame, metadata
