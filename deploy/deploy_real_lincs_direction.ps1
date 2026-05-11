param(
    [string]$Space = "Qiangaoqing/drugreflector-api"
)

$ErrorActionPreference = "Stop"

Write-Host "Uploading DrugReflector real LINCS direction backend to Hugging Face Space: $Space"

hf upload $Space drugreflector/api.py drugreflector/api.py --repo-type space --commit-message "Enable real L1000CDS2 direction evidence"
hf upload $Space drugreflector/directionality.py drugreflector/directionality.py --repo-type space --commit-message "Enable real L1000CDS2 direction evidence"

$env:DRUGREFLECTOR_HF_SPACE = $Space
@'
from huggingface_hub import HfApi
import os

space = os.environ.get("DRUGREFLECTOR_HF_SPACE", "Qiangaoqing/drugreflector-api")
HfApi().restart_space(space)
print(f"Restarted {space}")
'@ | python -

Write-Host "Deployment requested. Wait until the Space runtime is RUNNING, then test /api/health."
