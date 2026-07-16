#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
ba-arap-nma generate \
  --pdb input_state.pdb \
  --chain A \
  --modes 1 \
  --orientations 1 \
  --max-amplitude 1 \
  --amplitude-step 0.5 \
  --output output_generate \
  --overwrite
ba-arap-nma two-state \
  --input-pdb input_state.pdb \
  --target-pdb paired_second_state.pdb \
  --input-chain A \
  --target-chain A \
  --modes 1 \
  --orientations 1 \
  --analysis-modes 2 \
  --summary-k 1 \
  --min-identity 0 \
  --min-coverage 0 \
  --max-amplitude 1 \
  --amplitude-step 0.5 \
  --output output_two_state \
  --overwrite
printf '\nSmoke test completed successfully.\n'
