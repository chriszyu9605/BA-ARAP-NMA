$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot
ba-arap-nma two-state `
  --input-pdb input_state.pdb --target-pdb paired_second_state.pdb `
  --input-chain A --target-chain A `
  --modes 1 --orientations -1 `
  --analysis-modes 20 --summary-k 10 `
  --max-amplitude 5 --amplitude-step 0.25 `
  --output output `
  --overwrite
