"""Objective compound annotations from the Broad Repurposing Hub and PubChem."""

from __future__ import annotations

import io
import json
import os
import ssl
from functools import lru_cache
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd

APP_DIR = Path(__file__).resolve().parent
REPO_ROOT = APP_DIR.parent
REPO_CACHE_DIR = REPO_ROOT / ".repo_cache"
MAX_REPO_DOWNLOAD_BYTES = int(os.getenv("DRUGREFLECTOR_MAX_REPO_DOWNLOAD_MB", "100")) * 1024 * 1024

SAMPLE_INFO_URL = "https://repo-hub.broadinstitute.org/public/data/repo-sample-annotation-20240610.txt"
DRUG_INFO_URL = "https://repo-hub.broadinstitute.org/public/data/repo-drug-annotation-20200324.txt"
PUBCHEM_PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
    "{cid}/property/Title,IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES/JSON"
)
PUBCHEM_SYNONYMS_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"


def _canonical(value: object) -> str:
    return " ".join(str(value).strip().lower().split())


def _core_brd_id(value: object) -> str:
    text = str(value).strip().upper()
    parts = text.split("-")
    if len(parts) >= 3 and parts[0] == "BRD":
        return "-".join(parts[:2])
    return text


def _download_text(url: str, cache_name: str) -> str:
    REPO_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = REPO_CACHE_DIR / cache_name
    if cache_path.exists():
        return cache_path.read_text(encoding="utf-8", errors="replace")

    request = Request(url, headers={"User-Agent": "DrugReflector/1.0"})
    context = ssl.create_default_context()
    context.check_hostname = False
    context.verify_mode = ssl.CERT_NONE

    last_error: Exception | None = None
    for _ in range(3):
        try:
            with urlopen(request, timeout=120, context=context) as response:
                raw = response.read(MAX_REPO_DOWNLOAD_BYTES + 1)
            if len(raw) > MAX_REPO_DOWNLOAD_BYTES:
                limit_mb = MAX_REPO_DOWNLOAD_BYTES // (1024 * 1024)
                raise ValueError(f"Annotation file is too large to process safely: {url} (>{limit_mb} MB)")
            cache_path.write_bytes(raw)
            return raw.decode("utf-8", errors="replace")
        except (HTTPError, URLError) as exc:  # pragma: no cover - network dependent
            last_error = exc
            continue

    raise ValueError(f"Could not download Repurposing Hub annotation file: {url}") from last_error


def _read_repo_table(url: str, cache_name: str) -> pd.DataFrame:
    text = _download_text(url, cache_name)
    lines = [line for line in text.splitlines() if line and not line.startswith("!")]
    if not lines:
        raise ValueError(f"Repurposing Hub file {cache_name} did not contain a data table.")
    return pd.read_csv(io.StringIO("\n".join(lines)), sep="\t")


def _download_json(url: str, cache_name: str) -> dict[str, Any]:
    return json.loads(_download_text(url, cache_name))


def _prepare_sample_annotations() -> pd.DataFrame:
    sample = _read_repo_table(SAMPLE_INFO_URL, "repo-sample-annotation-20240610.txt")
    sample = sample.rename(columns={"broad_id": "full_broad_id"})
    sample["compound"] = sample["full_broad_id"].map(_core_brd_id)
    sample["pert_iname_norm"] = sample["pert_iname"].map(_canonical)
    sample["has_pubchem"] = sample["pubchem_cid"].notna().astype(int)
    sample["has_smiles"] = sample["smiles"].notna().astype(int)
    sample["purity_num"] = pd.to_numeric(sample["purity"], errors="coerce").fillna(-1.0)
    sample["deprecated_rank"] = sample["deprecated_broad_id"].fillna("").eq("").astype(int)
    sample = sample.sort_values(
        by=["compound", "deprecated_rank", "has_pubchem", "has_smiles", "purity_num"],
        ascending=[True, False, False, False, False],
    )
    sample = sample.drop_duplicates(subset=["compound"], keep="first")
    return sample


def _prepare_drug_annotations() -> pd.DataFrame:
    drug = _read_repo_table(DRUG_INFO_URL, "repo-drug-annotation-20200324.txt")
    drug["pert_iname_norm"] = drug["pert_iname"].map(_canonical)
    drug["annotation_score"] = drug[["clinical_phase", "moa", "target", "disease_area", "indication"]].notna().sum(axis=1)
    drug = drug.sort_values(by=["pert_iname_norm", "annotation_score"], ascending=[True, False])
    drug = drug.drop_duplicates(subset=["pert_iname_norm"], keep="first")
    return drug


@lru_cache(maxsize=1)
def _annotation_table() -> pd.DataFrame:
    sample = _prepare_sample_annotations()
    drug = _prepare_drug_annotations()
    merged = sample.merge(
        drug[["pert_iname_norm", "clinical_phase", "moa", "target", "disease_area", "indication"]],
        on="pert_iname_norm",
        how="left",
    )
    return merged


def _normalize_optional_text(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    return text or None


def _unique_synonyms(values: list[object], *, limit: int = 8) -> list[str]:
    seen: set[str] = set()
    synonyms: list[str] = []
    for value in values:
        text = _normalize_optional_text(value)
        if not text:
            continue
        key = text.casefold()
        if key in seen:
            continue
        seen.add(key)
        synonyms.append(text)
        if len(synonyms) >= limit:
            break
    return synonyms


@lru_cache(maxsize=4096)
def _pubchem_details(pubchem_cid: int) -> dict[str, Any]:
    try:
        prop_payload = _download_json(
            PUBCHEM_PROPERTY_URL.format(cid=pubchem_cid),
            f"pubchem-properties-{pubchem_cid}.json",
        )
    except Exception:  # pragma: no cover - network dependent
        return {}

    properties = prop_payload.get("PropertyTable", {}).get("Properties", [])
    prop = properties[0] if properties else {}

    try:
        synonym_payload = _download_json(
            PUBCHEM_SYNONYMS_URL.format(cid=pubchem_cid),
            f"pubchem-synonyms-{pubchem_cid}.json",
        )
        synonym_records = synonym_payload.get("InformationList", {}).get("Information", [])
        raw_synonyms = synonym_records[0].get("Synonym", []) if synonym_records else []
    except Exception:  # pragma: no cover - network dependent
        raw_synonyms = []

    return {
        "pubchem_title": _normalize_optional_text(prop.get("Title")),
        "chemical_name": _normalize_optional_text(prop.get("IUPACName")),
        "molecular_formula": _normalize_optional_text(prop.get("MolecularFormula")),
        "molecular_weight": (
            float(prop["MolecularWeight"])
            if prop.get("MolecularWeight") is not None
            else None
        ),
        "canonical_smiles": _normalize_optional_text(prop.get("CanonicalSMILES")),
        "isomeric_smiles": _normalize_optional_text(prop.get("IsomericSMILES")),
        "synonyms": _unique_synonyms(raw_synonyms),
    }


def get_compound_annotations(compounds: list[str]) -> dict[str, dict[str, Any]]:
    if not compounds:
        return {}

    table = _annotation_table()
    wanted = {_core_brd_id(compound) for compound in compounds}
    subset = table[table["compound"].isin(wanted)].copy()
    annotations: dict[str, dict[str, Any]] = {}

    for _, row in subset.iterrows():
        cid = row.get("pubchem_cid")
        pubchem_cid = int(cid) if pd.notna(cid) else None
        pubchem = _pubchem_details(pubchem_cid) if pubchem_cid is not None else {}
        annotations[str(row["compound"])] = {
            "display_name": _normalize_optional_text(row.get("pert_iname")),
            "full_broad_id": _normalize_optional_text(row.get("full_broad_id")),
            "clinical_phase": _normalize_optional_text(row.get("clinical_phase")),
            "moa": _normalize_optional_text(row.get("moa")),
            "target": _normalize_optional_text(row.get("target")),
            "disease_area": _normalize_optional_text(row.get("disease_area")),
            "indication": _normalize_optional_text(row.get("indication")),
            "vendor_name": _normalize_optional_text(row.get("vendor_name")),
            "smiles": _normalize_optional_text(row.get("smiles")),
            "inchikey": _normalize_optional_text(row.get("InChIKey")),
            "trade_name": None,
            "pubchem_title": _normalize_optional_text(pubchem.get("pubchem_title")),
            "chemical_name": _normalize_optional_text(pubchem.get("chemical_name")),
            "molecular_formula": _normalize_optional_text(pubchem.get("molecular_formula")),
            "molecular_weight": pubchem.get("molecular_weight"),
            "canonical_smiles": _normalize_optional_text(pubchem.get("canonical_smiles")),
            "isomeric_smiles": _normalize_optional_text(pubchem.get("isomeric_smiles")),
            "synonyms": pubchem.get("synonyms", []),
            "pubchem_cid": pubchem_cid,
            "pubchem_url": (
                f"https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem_cid}"
                if pubchem_cid is not None
                else None
            ),
            "structure_image": (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/PNG?image_size=large"
                if pubchem_cid is not None
                else None
            ),
            "source": "Broad Repurposing Hub + PubChem",
        }

    return annotations
