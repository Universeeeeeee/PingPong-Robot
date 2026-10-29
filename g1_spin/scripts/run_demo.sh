#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
.venv/bin/python -m pytest -q
.venv/bin/python -m g1_spin.cli synthetic-benchmark --output "${1:-outputs/demo}"
