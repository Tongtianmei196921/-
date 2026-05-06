$ErrorActionPreference = "Stop"

$root = Resolve-Path (Join-Path $PSScriptRoot "..\..\..")
Set-Location $root
New-Item -ItemType Directory -Force -Path "checkpoints" | Out-Null

$required = @(
  "checkpoints\model_fold_0.pt",
  "checkpoints\model_fold_1.pt",
  "checkpoints\model_fold_2.pt"
)

if (($required | Where-Object { -not (Test-Path $_) }).Count -eq 0) {
  Write-Host "Checkpoint files already exist."
  exit 0
}

python -m pip install --upgrade pip zenodo-get
zenodo_get --output-dir checkpoints 16912444
