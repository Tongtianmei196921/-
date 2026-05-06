#!/usr/bin/env sh
set -eu

cd "$(dirname "$0")/../../.."
mkdir -p checkpoints

if [ -f checkpoints/model_fold_0.pt ] && [ -f checkpoints/model_fold_1.pt ] && [ -f checkpoints/model_fold_2.pt ]; then
  echo "Checkpoint files already exist."
  exit 0
fi

python -m pip install --upgrade pip zenodo-get
zenodo_get --output-dir checkpoints 16912444
