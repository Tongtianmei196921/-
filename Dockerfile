FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

COPY pyproject.toml README.md LICENSE.txt ./
COPY drugreflector ./drugreflector
COPY signature_refinement ./signature_refinement
COPY checkpoints ./checkpoints

RUN pip install --upgrade pip && \
    pip install ".[api]"

EXPOSE 8000

CMD ["python", "-m", "uvicorn", "drugreflector.api:app", "--host", "0.0.0.0", "--port", "8000"]
