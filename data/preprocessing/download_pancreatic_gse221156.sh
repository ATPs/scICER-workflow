#!/usr/bin/env bash
set -euo pipefail

out_dir="${1:-data/raw/pancreas_gse221156}"
url="https://datasets.cellxgene.cziscience.com/74a82e6d-2c1a-4328-9011-c6a437747854.h5ad"
file_name="74a82e6d-2c1a-4328-9011-c6a437747854.h5ad"

mkdir -p "$out_dir"
cd "$out_dir"

if command -v curl >/dev/null 2>&1; then
  curl -L --continue-at - -o "$file_name" "$url"
else
  wget -c -O "$file_name" "$url"
fi
