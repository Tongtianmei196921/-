"""FastAPI service for DrugReflector inference and UI."""

from __future__ import annotations

import io
import os
import tempfile
from pathlib import Path
from threading import Lock

import numpy as np
import pandas as pd
from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from .directionality import get_direction_evidence
from .drug_reflector import DrugReflector
from .geo import build_geo_signature, preview_geo
from .prepare import prepare_uploaded_input
from .report import build_docx_report
from .repo_annotations import get_compound_annotations
from .utils import load_h5ad_file

APP_DIR = Path(__file__).resolve().parent
REPO_ROOT = APP_DIR.parent
STATIC_DIR = APP_DIR / "static"
FRONTEND_DIST_DIR = REPO_ROOT / "frontend" / "dist"
FRONTEND_ASSETS_DIR = FRONTEND_DIST_DIR / "assets"
CHECKPOINT_NAMES = ["model_fold_0.pt", "model_fold_1.pt", "model_fold_2.pt"]
DEFAULT_CORS_ORIGINS = [
    "https://drugreflector.vercel.app",
    "http://localhost:5173",
    "http://127.0.0.1:5173",
    "http://localhost:8000",
    "http://127.0.0.1:8000",
]
MAX_UPLOAD_BYTES = int(os.getenv("DRUGREFLECTOR_MAX_UPLOAD_MB", "200")) * 1024 * 1024
MAX_REPORT_SAMPLES = int(os.getenv("DRUGREFLECTOR_MAX_REPORT_SAMPLES", "5"))
MAX_REPORT_ROWS = int(os.getenv("DRUGREFLECTOR_MAX_REPORT_ROWS", "1000"))
MIN_N_TOP = 1
MAX_N_TOP = 500


def _default_checkpoint_paths() -> list[str]:
    base_dir = REPO_ROOT / "checkpoints"
    return [str(base_dir / name) for name in CHECKPOINT_NAMES]


def resolve_checkpoint_paths() -> list[str]:
    raw_paths = os.getenv("DRUGREFLECTOR_MODEL_PATHS", "").strip()
    if raw_paths:
        separators = [";", ","]
        parts = [raw_paths]
        for separator in separators:
            if separator in raw_paths:
                parts = [segment.strip() for segment in raw_paths.split(separator)]
                break
        paths = [segment for segment in parts if segment]
    else:
        checkpoint_dir = os.getenv("DRUGREFLECTOR_CHECKPOINT_DIR", "").strip()
        if checkpoint_dir:
            base_dir = Path(checkpoint_dir)
            paths = [str(base_dir / name) for name in CHECKPOINT_NAMES]
        else:
            paths = _default_checkpoint_paths()

    if len(paths) != 3:
        raise RuntimeError(
            "DrugReflector requires exactly 3 checkpoint paths. "
            "Set DRUGREFLECTOR_MODEL_PATHS or DRUGREFLECTOR_CHECKPOINT_DIR."
        )
    return paths


def resolve_checkpoint_dir() -> Path:
    paths = [Path(path).resolve() for path in resolve_checkpoint_paths()]
    parents = {path.parent for path in paths}
    if len(parents) == 1:
        return next(iter(parents))
    return REPO_ROOT / "checkpoints"


def checkpoint_status() -> dict[str, object]:
    paths = [Path(path).resolve() for path in resolve_checkpoint_paths()]
    status = {
        CHECKPOINT_NAMES[idx]: paths[idx].exists()
        for idx in range(min(len(paths), len(CHECKPOINT_NAMES)))
    }
    return {
        "checkpoints": status,
        "all_ready": all(status.values()),
        "checkpoint_dir": str(resolve_checkpoint_dir()),
    }


class SignatureInput(BaseModel):
    name: str = Field(default="signature")
    scores: dict[str, float] = Field(min_length=1)


class PredictRequest(BaseModel):
    signatures: list[SignatureInput] = Field(min_length=1)
    n_top: int = Field(default=50, ge=1, le=500)


class GeoPredictRequest(BaseModel):
    accession: str = Field(min_length=4)
    n_top: int = Field(default=50, ge=1, le=500)
    group_column: str | None = None
    group1_value: str | None = None
    group2_value: str | None = None
    control_keyword: str | None = None
    case_keyword: str | None = None


class ReportRequest(BaseModel):
    response: dict[str, object]
    sample: str | None = None
    locale: str = Field(default="zh", max_length=8)


_MODEL_LOCK = Lock()
_MODEL: DrugReflector | None = None
_MODEL_ERROR: str | None = None


def get_model() -> DrugReflector:
    global _MODEL, _MODEL_ERROR

    if _MODEL is not None:
        return _MODEL

    with _MODEL_LOCK:
        if _MODEL is not None:
            return _MODEL

        try:
            checkpoint_paths = resolve_checkpoint_paths()
            missing_paths = [path for path in checkpoint_paths if not Path(path).exists()]
            if missing_paths:
                raise FileNotFoundError(
                    "Missing checkpoint files: " + ", ".join(missing_paths)
                )

            _MODEL = DrugReflector(checkpoint_paths=checkpoint_paths)
            _MODEL_ERROR = None
            return _MODEL
        except Exception as exc:  # pragma: no cover
            _MODEL_ERROR = str(exc)
            raise


def _cors_origins() -> list[str]:
    raw_origins = os.getenv("DRUGREFLECTOR_CORS_ORIGINS", "").strip()
    if not raw_origins:
        return DEFAULT_CORS_ORIGINS
    return [origin.strip() for origin in raw_origins.split(",") if origin.strip()]


def _validate_n_top(n_top: int) -> int:
    try:
        value = int(n_top)
    except (TypeError, ValueError) as exc:
        raise HTTPException(status_code=422, detail="n_top must be an integer.") from exc
    if value < MIN_N_TOP or value > MAX_N_TOP:
        raise HTTPException(
            status_code=422,
            detail=f"n_top must be between {MIN_N_TOP} and {MAX_N_TOP}.",
        )
    return value


async def _read_upload_limited(file: UploadFile) -> bytes:
    content = await file.read(MAX_UPLOAD_BYTES + 1)
    if len(content) > MAX_UPLOAD_BYTES:
        limit_mb = MAX_UPLOAD_BYTES // (1024 * 1024)
        raise HTTPException(
            status_code=413,
            detail=f"Uploaded file is too large. Maximum allowed size is {limit_mb} MB.",
        )
    return content


def _validate_report_request(request: ReportRequest) -> str:
    locale = request.locale.lower().strip()
    if locale not in {"zh", "en"}:
        raise HTTPException(status_code=422, detail="locale must be 'zh' or 'en'.")

    data = request.response.get("data")
    if not isinstance(data, dict):
        raise HTTPException(status_code=422, detail="Report response must contain data.")

    samples = data.get("samples") or []
    if not isinstance(samples, list) or len(samples) > MAX_REPORT_SAMPLES:
        raise HTTPException(
            status_code=413,
            detail=f"Report request is too large. Maximum samples: {MAX_REPORT_SAMPLES}.",
        )

    results = data.get("results") or {}
    if not isinstance(results, dict):
        raise HTTPException(status_code=422, detail="Report response results are invalid.")

    row_count = 0
    for rows in results.values():
        if not isinstance(rows, list):
            raise HTTPException(status_code=422, detail="Report result rows are invalid.")
        row_count += len(rows)
    if row_count > MAX_REPORT_ROWS:
        raise HTTPException(
            status_code=413,
            detail=f"Report request is too large. Maximum result rows: {MAX_REPORT_ROWS}.",
        )

    return locale


def _predictions_to_payload(
    predictions: pd.DataFrame,
    annotations: dict[str, dict[str, object]] | None = None,
    direction_evidence: dict[str, dict[str, dict[str, object]]] | None = None,
) -> dict[str, list[dict[str, float | int | str]]]:
    payload: dict[str, list[dict[str, float | int | str]]] = {}
    annotations = annotations or {}
    direction_evidence = direction_evidence or {}

    signature_names = list(dict.fromkeys(column[1] for column in predictions.columns))
    for signature_name in signature_names:
        ordered = predictions[("prob", signature_name)].sort_values(ascending=False)
        rows = []
        for compound in ordered.index:
            rows.append(
                {
                    "compound": compound,
                    "rank": int(predictions.loc[compound, ("rank", signature_name)]) + 1,
                    "logit": float(predictions.loc[compound, ("logit", signature_name)]),
                    "prob": float(predictions.loc[compound, ("prob", signature_name)]),
                    "annotation": annotations.get(str(compound)),
                    "direction": direction_evidence.get(signature_name, {}).get(str(compound)),
                }
            )
        payload[signature_name] = rows

    return payload


def _frame_to_adata(frame: pd.DataFrame):
    import anndata

    return anndata.AnnData(
        X=frame.values.astype(np.float32),
        obs=pd.DataFrame(index=frame.index),
        var=pd.DataFrame(index=frame.columns),
    )


def _sample_gene_names() -> list[str]:
    try:
        model = get_model()
        return list(model.model.dimensions["var_names"][0])
    except Exception:
        return [f"GENE_{idx:04d}" for idx in range(978)]


def _normalize_gene_name(value: object) -> str:
    text = str(value).strip().upper()
    if not text:
        return ""
    text = text.split(".")[0]
    if text.endswith("_AT"):
        text = text[:-3]
    return "".join(ch for ch in text if ch.isalnum() or ch == "-")


def _input_quality_summary(data) -> dict[str, float | int]:
    if hasattr(data, "var_names"):
        input_genes = [_normalize_gene_name(gene) for gene in data.var_names]
    else:
        input_genes = [_normalize_gene_name(gene) for gene in data.columns]

    input_gene_set = {gene for gene in input_genes if gene}
    model_gene_set = {_normalize_gene_name(gene) for gene in _sample_gene_names()}
    model_gene_set.discard("")

    overlap = input_gene_set & model_gene_set
    missing = model_gene_set - input_gene_set
    ratio = (len(overlap) / len(model_gene_set)) if model_gene_set else 0.0

    return {
        "input_gene_count": int(len(input_gene_set)),
        "model_gene_count": int(len(model_gene_set)),
        "overlap_gene_count": int(len(overlap)),
        "missing_model_gene_count": int(len(missing)),
        "overlap_ratio": float(ratio),
    }


def _generate_sample_csv() -> str:
    rng = np.random.default_rng(42)
    genes = _sample_gene_names()
    df = pd.DataFrame(
        rng.standard_normal((3, len(genes))).astype(np.float32),
        index=["Sample_A", "Sample_B", "Sample_C"],
        columns=genes,
    )
    buffer = io.StringIO()
    df.to_csv(buffer)
    return buffer.getvalue()


async def _uploaded_file_to_adata(file: UploadFile):
    filename = file.filename or ""
    extension = Path(filename).suffix.lower()
    content = await _read_upload_limited(file)

    try:
        if extension == ".h5ad":
            suffix = extension or ".h5ad"
            with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as temp_file:
                temp_path = Path(temp_file.name)
                temp_file.write(content)
            try:
                adata = load_h5ad_file(str(temp_path))
            finally:
                temp_path.unlink(missing_ok=True)
            return adata

        if extension in {".csv", ".tsv"}:
            separator = "\t" if extension == ".tsv" else ","
            frame = pd.read_csv(io.BytesIO(content), sep=separator, index_col=0)
            return _frame_to_adata(frame)

    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=422, detail=f"Could not parse file: {exc}") from exc

    raise HTTPException(status_code=400, detail=f"Unsupported file type: {extension!r}")


async def _prepare_uploaded_file(
    file: UploadFile,
    *,
    group_column: str | None = None,
    group1_value: str | None = None,
    group2_value: str | None = None,
    sample_id_column: str | None = None,
):
    filename = file.filename or ""
    content = await _read_upload_limited(file)
    if not filename:
        raise HTTPException(status_code=400, detail="Uploaded file must have a filename.")

    try:
        return prepare_uploaded_input(
            filename,
            content,
            group_column=group_column,
            group1_value=group1_value,
            group2_value=group2_value,
            sample_id_column=sample_id_column,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


def _friendly_geo_error(exc: ValueError, accession: str) -> str:
    message = str(exc)
    accession = accession.upper()
    if "unknown url type" in message and "NONE" in message:
        return (
            f"{accession} 的 GEO 记录里包含 NONE 这类空的补充文件链接，"
            "系统已判定它不是可下载的表达矩阵。请更换 GEO 编号，"
            "或先把该数据集整理成 h5ad / CSV 表达矩阵后再上传。"
        )
    if (
        "does not expose an expression matrix" in message
        or "no per-sample supplementary expression files" in message
    ):
        if (
            "single-cell viewer/archive" in message
            or ".cloupe" in message
            or "10x analysis archives" in message
        ):
            return (
                f"{accession} 暂时不能直接自动分析：该 GEO 没有提供可直接使用的基因表达矩阵，"
                "而是提供了 .cloupe、.tar、.h5 或 10x analysis archive 这类单细胞原始/浏览器文件。"
                "这些文件需要先经过单细胞流程整理成 h5ad、伪 bulk 表达矩阵或基因级差异表后再上传。"
            )
        return (
            f"{accession} 暂时不能自动解析：series matrix 中没有可用表达矩阵，"
            "也没有逐样本表达补充文件。请换一个包含表达矩阵的 GEO，"
            "或下载原始/处理后数据整理成 h5ad / CSV 后上传。"
        )
    if "differential miRNA table" in message or "miRNA IDs cannot be objectively treated as genes" in message:
        return (
            f"{accession} 找到了作者提供的差异 miRNA 表，但这不是普通基因差异表达表。"
            "DrugReflector 需要人类基因符号作为输入，不能把 miRNA ID 直接当作基因。"
            "如需分析该数据，请先基于可靠数据库完成 miRNA 靶基因映射并形成基因级 signature，"
            "或上传整理好的基因差异表。"
        )
    return message


def _public_geo_error(exc: ValueError, accession: str) -> str:
    message = str(exc)
    accession = accession.upper()
    if "unknown url type" in message and "NONE" in message:
        return (
            f"{accession} 的 GEO 记录里包含 NONE 这类空的补充文件链接，"
            "系统已判定它不是可下载的表达矩阵。请更换 GEO 编号，"
            "或先把该数据集整理成 h5ad / CSV 表达矩阵后再上传。"
        )
    if (
        "does not expose an expression matrix" in message
        or "no per-sample supplementary expression files" in message
    ):
        return (
            f"{accession} 暂时不能自动解析：series matrix 中没有可用表达矩阵，"
            "也没有逐样本表达补充文件。请换一个包含表达矩阵的 GEO，"
            "或下载原始 / 处理后数据整理成 h5ad / CSV 后上传。"
        )
    if "differential miRNA table" in message or "miRNA IDs cannot be objectively treated as genes" in message:
        return (
            f"{accession} 找到了作者提供的差异 miRNA 表，但这不是普通基因差异表达表。"
            "DrugReflector 需要人类基因符号作为输入，不能把 miRNA ID 直接当作基因。"
            "如需分析该数据，请先基于可靠数据库完成 miRNA 靶基因映射并形成基因级 signature，"
            "或上传整理好的基因差异表。"
        )
    return message


app = FastAPI(
    title="DrugReflector API",
    version="1.0.0",
    description="Inference service for compound ranking from gene expression signatures.",
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins(),
    allow_methods=["*"],
    allow_headers=["*"],
)
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
if FRONTEND_ASSETS_DIR.exists():
    app.mount("/assets", StaticFiles(directory=FRONTEND_ASSETS_DIR), name="assets")


@app.get("/", include_in_schema=False)
def index() -> FileResponse:
    if (FRONTEND_DIST_DIR / "index.html").exists():
        return FileResponse(FRONTEND_DIST_DIR / "index.html")
    return FileResponse(STATIC_DIR / "index.html")


@app.get("/health")
def health() -> dict[str, object]:
    response = {"status": "ok", **checkpoint_status()}
    try:
        model = get_model()
        response["n_compounds"] = model.n_compounds
        response["checkpoints"] = model.checkpoint_paths
    except Exception:
        response["status"] = "error"
        response["detail"] = _MODEL_ERROR or "Model failed to load."
    return response


@app.get("/metadata")
def metadata() -> dict[str, object]:
    model = get_model()
    input_gene_counts = [len(genes) for genes in model.model.dimensions["var_names"]]
    return {
        "n_compounds": model.n_compounds,
        "checkpoint_paths": model.checkpoint_paths,
        "input_gene_counts": input_gene_counts,
        "output_preview": model.compound_names[:10],
    }


@app.get("/api/health")
def api_health() -> dict[str, object]:
    return {"status": "ok", **checkpoint_status()}


@app.get("/api/checkpoints")
def api_checkpoints() -> dict[str, object]:
    return checkpoint_status()


@app.get("/api/sample")
def api_sample() -> StreamingResponse:
    csv_text = _generate_sample_csv()
    return StreamingResponse(
        io.BytesIO(csv_text.encode()),
        media_type="text/csv",
        headers={"Content-Disposition": "attachment; filename=sample_vscores.csv"},
    )


@app.post("/api/prepare-input")
async def api_prepare_input(
    file: UploadFile = File(...),
    group_column: str | None = Form(default=None),
    group1_value: str | None = Form(default=None),
    group2_value: str | None = Form(default=None),
    sample_id_column: str | None = Form(default=None),
) -> dict[str, object]:
    prepared = await _prepare_uploaded_file(
        file,
        group_column=group_column,
        group1_value=group1_value,
        group2_value=group2_value,
        sample_id_column=sample_id_column,
    )
    return {"success": True, "preparation": prepared.summary}


@app.get("/api/geo/preview")
def api_geo_preview(accession: str) -> dict[str, object]:
    try:
        return preview_geo(accession)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=_public_geo_error(exc, accession)) from exc


@app.post("/api/geo/predict")
def api_geo_predict(request: GeoPredictRequest) -> dict[str, object]:
    status = checkpoint_status()
    if not status["all_ready"]:
        checkpoint_dir = status["checkpoint_dir"]
        raise HTTPException(
            status_code=503,
            detail=(
                "Model checkpoints not found. "
                f'Download via: zenodo_get --output-dir "{checkpoint_dir}" 16912444'
            ),
        )

    try:
        frame, geo_meta = build_geo_signature(
            request.accession,
            group_column=request.group_column,
            group1_value=request.group1_value,
            group2_value=request.group2_value,
            control_keyword=request.control_keyword,
            case_keyword=request.case_keyword,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=_public_geo_error(exc, request.accession)) from exc

    model = get_model()
    predictions = model.predict(frame, n_top=request.n_top)
    annotations = get_compound_annotations(predictions.index.tolist())
    compounds_by_sample = {
        str(sample_name): list(predictions[("prob", str(sample_name))].sort_values(ascending=False).index)
        for sample_name in frame.index
    }
    direction_evidence = get_direction_evidence(frame, compounds_by_sample)
    input_quality = _input_quality_summary(frame)
    payload = _predictions_to_payload(
        predictions,
        annotations=annotations,
        direction_evidence=direction_evidence,
    )
    sample_names = list(payload.keys())

    return {
        "success": True,
        "data": {
            "samples": sample_names,
            "results": payload,
        },
        "meta": {
            "source": "geo",
            "filename": request.accession.upper(),
            "n_samples": 1,
            "n_compounds": int(len(predictions)),
            "n_genes": int(frame.shape[1]),
            "input_quality": input_quality,
            **geo_meta,
        },
    }


@app.post("/api/report/docx")
def api_report_docx(request: ReportRequest) -> StreamingResponse:
    locale = _validate_report_request(request)
    try:
        content = build_docx_report(
            request.response,
            sample=request.sample,
            locale=locale,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail="Could not generate report.") from exc

    sample_name = (request.sample or "report").replace("/", "_").replace("\\", "_")
    filename = f"drugreflector_{sample_name}_report.docx"
    return StreamingResponse(
        io.BytesIO(content),
        media_type="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


@app.post("/predict/vscores")
def predict_vscores(request: PredictRequest) -> dict[str, object]:
    model = get_model()

    frame = pd.DataFrame(
        [signature.scores for signature in request.signatures],
        index=[signature.name for signature in request.signatures],
    ).fillna(0.0)

    predictions = model.predict(frame, n_top=request.n_top)
    annotations = get_compound_annotations(predictions.index.tolist())
    compounds_by_sample = {
        str(sample_name): list(predictions[("prob", str(sample_name))].sort_values(ascending=False).index)
        for sample_name in frame.index
    }
    direction_evidence = get_direction_evidence(frame, compounds_by_sample)
    input_quality = _input_quality_summary(frame)
    return {
        "n_signatures": len(request.signatures),
        "n_top": request.n_top,
        "input_quality": input_quality,
        "results": _predictions_to_payload(
            predictions,
            annotations=annotations,
            direction_evidence=direction_evidence,
        ),
    }


@app.post("/predict/h5ad")
async def predict_h5ad(
    file: UploadFile = File(...),
    n_top: int = Form(default=50),
) -> dict[str, object]:
    if not file.filename:
        raise HTTPException(status_code=400, detail="Uploaded file must have a filename.")

    n_top = _validate_n_top(n_top)
    adata = await _uploaded_file_to_adata(file)
    model = get_model()
    predictions = model.predict(adata, n_top=n_top)
    annotations = get_compound_annotations(predictions.index.tolist())
    compounds_by_sample = {
        str(sample_name): list(predictions[("prob", str(sample_name))].sort_values(ascending=False).index)
        for sample_name in adata.obs_names.astype(str)
    }
    direction_evidence = get_direction_evidence(adata, compounds_by_sample)
    input_quality = _input_quality_summary(adata)
    return {
        "input_file": file.filename,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_top": n_top,
        "input_quality": input_quality,
        "results": _predictions_to_payload(
            predictions,
            annotations=annotations,
            direction_evidence=direction_evidence,
        ),
    }


@app.post("/api/predict")
async def api_predict(
    file: UploadFile = File(...),
    n_top: int = Form(default=50),
    group_column: str | None = Form(default=None),
    group1_value: str | None = Form(default=None),
    group2_value: str | None = Form(default=None),
    sample_id_column: str | None = Form(default=None),
) -> dict[str, object]:
    n_top = _validate_n_top(n_top)
    status = checkpoint_status()
    if not status["all_ready"]:
        checkpoint_dir = status["checkpoint_dir"]
        raise HTTPException(
            status_code=503,
            detail=(
                "Model checkpoints not found. "
                f'Download via: zenodo_get --output-dir "{checkpoint_dir}" 16912444'
            ),
        )

    prepared = await _prepare_uploaded_file(
        file,
        group_column=group_column,
        group1_value=group1_value,
        group2_value=group2_value,
        sample_id_column=sample_id_column,
    )
    if prepared.adata is None:
        raise HTTPException(
            status_code=422,
            detail="Input could not be prepared objectively. Review the preparation summary first.",
        )

    adata = prepared.adata
    model = get_model()
    predictions = model.predict(adata, n_top=n_top)
    annotations = get_compound_annotations(predictions.index.tolist())
    compounds_by_sample = {
        str(sample_name): list(predictions[("prob", str(sample_name))].sort_values(ascending=False).index)
        for sample_name in adata.obs_names.astype(str)
    }
    direction_evidence = get_direction_evidence(adata, compounds_by_sample)
    input_quality = _input_quality_summary(adata)
    payload = _predictions_to_payload(
        predictions,
        annotations=annotations,
        direction_evidence=direction_evidence,
    )
    sample_names = list(payload.keys())

    return {
        "success": True,
        "data": {
            "samples": sample_names,
            "results": payload,
        },
        "meta": {
            "n_samples": int(adata.n_obs),
            "n_compounds": int(len(predictions)),
            "n_genes": int(adata.n_vars),
            "filename": file.filename,
            "input_quality": input_quality,
            "preparation": prepared.summary,
        },
    }
