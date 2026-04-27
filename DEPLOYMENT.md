# DrugReflector Deployment

This project is easiest to deploy as a lightweight HTTP inference service for internal research use.

If you just want the local research UI on Windows, the fastest option is:

```powershell
.\start_ui.ps1
```

or double-click:

```text
start_ui.bat
```

The launcher prefers port `8000`, but if that port is already occupied it will automatically choose another local port such as `8010` and print the final URL.

If you are editing the React UI in `frontend/`, build it with:

```bash
cd frontend
npm install
npm run build
```

## Option 1: Run directly with Python

Install the package and API dependencies:

```bash
pip install -e ".[api]"
```

Start the service:

```bash
drugreflector serve --host 0.0.0.0 --port 8000
```

Open the UI:

```bash
http://localhost:8000/
```

If your checkpoints live somewhere else, set one of these before starting:

```bash
export DRUGREFLECTOR_CHECKPOINT_DIR=/path/to/checkpoints
```

or

```bash
export DRUGREFLECTOR_MODEL_PATHS=/path/model_fold_0.pt,/path/model_fold_1.pt,/path/model_fold_2.pt
```

## Option 2: Run with Docker

Build the image:

```bash
docker build -t drugreflector-api .
```

Run the container:

```bash
docker run --rm -p 8000:8000 drugreflector-api
```

If you want to mount external checkpoints:

```bash
docker run --rm -p 8000:8000 \
  -e DRUGREFLECTOR_CHECKPOINT_DIR=/models \
  -v /your/checkpoints:/models \
  drugreflector-api
```

## Endpoints

- `GET /`: browser UI for researchers
- `GET /health`: health check and checkpoint status
- `GET /metadata`: model metadata
- `GET /api/checkpoints`: checkpoint readiness for the React UI
- `GET /api/sample`: sample CSV generator
- `POST /api/predict`: upload-based prediction endpoint used by the React UI
- `POST /predict/vscores`: submit JSON v-score signatures
- `POST /predict/h5ad`: upload an `.h5ad` file

## Example request

```bash
curl -X POST http://localhost:8000/predict/vscores \
  -H "Content-Type: application/json" \
  -d '{
    "n_top": 5,
    "signatures": [
      {
        "name": "demo",
        "scores": {
          "TP53": 1.2,
          "EGFR": -0.4,
          "CDKN1A": 0.8
        }
      }
    ]
  }'
```

Research UI is available at `http://localhost:8000/`.

Interactive docs are available at `http://localhost:8000/docs`.

There is also a Chinese quick-start guide in `QUICKSTART.zh-CN.md`.
