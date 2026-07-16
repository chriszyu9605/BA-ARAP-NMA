# BA-ARAP-NMA 1.0.0

**BA-ARAP-NMA** (backbone-angle-compatible as-rigid-as-possible normal-mode analysis) generates amplitude-ordered protein structures from low-frequency C-alpha elastic-network modes while reducing local C-alpha backbone distortion.

The release provides two user workflows:

- `generate`: target-free structure generation from one PDB file.
- `two-state`: structure generation plus retrospective comparison with a paired second-state PDB and modal-space analyses.

The paired second state is never used to generate BA-ARAP-NMA coordinates. It is used only for residue mapping, rigid alignment, cumulative-overlap analysis, displacement retention/concentration analysis, RMSD, transition coverage, and selection of the closest generated structures.

## Requirements

- Python 3.10 or later
- NumPy and SciPy, installed automatically
- PyMOL is optional and is used only to open generated `.pml` visualization scripts

## Installation

```bash
python -m venv .venv
source .venv/bin/activate       # Windows PowerShell: .\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install .
ba-arap-nma --version
```

For the exact Python 3.13.5 environment used to validate this release:

```bash
python -m pip install -r requirements-lock.txt
python -m pip install --no-deps .
python -m pip install "pytest>=8"
pytest -q
```

For Python 3.10–3.12, use the standard installation above; `pyproject.toml` resolves compatible NumPy and SciPy versions for the selected Python version.

## Quick smoke test

After installation:

```bash
cd examples/smoke_test
bash run_smoke_test.sh
```

Windows PowerShell:

```powershell
cd examples/smoke_test
.\run_smoke_test.ps1
```

## Input

BA-ARAP-NMA accepts PDB files. Use `--chain` for a single-state calculation, or `--input-chain` and `--target-chain` for a two-state calculation. Protein C-alpha atoms are extracted from the selected chain. Standard protein `ATOM` records and the protein-like modified residues MSE, SEC, and PYL are accepted.

## Target-free structure generation

```bash
ba-arap-nma generate \
  --pdb input.pdb \
  --chain A \
  --modes 1-10 \
  --orientations -1 1 \
  --max-amplitude 20 \
  --amplitude-step 0.25 \
  --output output_single_state
```

For each mode rank and orientation, the initial constrained mode is calculated once and held fixed across all requested amplitudes. With the default 0–20 Å range and 0.25 Å interval, each branch contains 81 frames including the input structure at frame 0.

The program refuses to write into a non-empty output directory. Choose a new directory or add `--overwrite` when replacement is intentional. Use `--quiet` to suppress progress messages.

## Two-state analysis

```bash
ba-arap-nma two-state \
  --input-pdb input_state.pdb \
  --target-pdb paired_second_state.pdb \
  --input-chain A \
  --target-chain A \
  --modes 1-10 \
  --orientations -1 1 \
  --analysis-modes 50 \
  --summary-k 10 \
  --max-amplitude 20 \
  --amplitude-step 0.25 \
  --output output_two_state
```

This workflow generates matched BA-ARAP-NMA and Linear C-alpha ANM branches. It also writes:

- cumulative overlap for every `K` from 1 to `--analysis-modes`;
- complete-space displacement retention (`eta`) and maximum attainable overlap (`sqrt(eta)`);
- the fraction of retained displacement carried by the first `K` modes (`rho_K`);
- the first-order local-distance change in the observed two-state displacement (`s_loc`);
- minimum-RMSD structures, transition coverage, local-angle error, PDB files, and a PyMOL comparison script.

## HCV helicase example

The included article-scale example uses the HCV helicase pair, which showed the clearest minimum-RMSD improvement in the associated study. The PDB files contain the common mapped C-alpha residues used for the two-state comparison.

```bash
cd examples/hcv_helicase
bash run_example.sh
```

Windows PowerShell:

```powershell
cd examples/hcv_helicase
.\run_example.ps1
```

The shortened example samples mode rank 1, orientation -1, and amplitudes 0–5 Å. It is a real protein example and can require several minutes or longer and substantial memory, depending on protein size, CPU, and the installed linear-algebra backend. Run the smoke test first to verify installation. Archived expected outputs are included in `expected_results/`.

## Output files

Each generated branch contains:

- `ca_structures.pdb`: amplitude-ordered multi-model C-alpha structures;
- `ca_trace_structures.pdb`: the same structures with sequential `CONECT` records;
- `view_ca_structures.pml`: PyMOL visualization script;
- `selected_structures/`: structures at the beginning, quarter points, midpoint, three-quarter point, and end;
- `frame_metrics.tsv`: the local-geometry measurements used in the study, plus target RMSD and transition coverage when a target is supplied;
- `branch_summary.json` and `run_config.json`;
- `fullatom_decorated_structures.pdb`, `selected_fullatom_structures/`, and `view_fullatom_decorated.pml` when the input PDB contains non-C-alpha protein atoms.

The decorated full-atom set moves the input atom records according to local C-alpha changes. It is intended for visualization and as a starting point for rebuilding and refinement; it is not an energy-minimized all-atom model.

At the run level, `branch_summary.tsv` lists all branches. The two-state workflow additionally writes `modal_space_analysis.tsv`, `modal_space_summary.tsv`, `best_structures/`, and `best_structure_comparison.tsv`.

## Metrics

`frame_metrics.tsv` reports only the measurements used in the associated study:

- harmonic virtual-bond distortion proxy (`k = 500`);
- adjacent C-alpha distance mean absolute error;
- adjacent-distance outlier counts above 0.5 Å and 1.0 Å;
- virtual-angle mean absolute error;
- virtual-angle outlier count above 10°;
- aligned target RMSD and transition coverage in the two-state workflow.

## Reproducibility and scope

BA-ARAP-NMA is a C-alpha structure-generation method, not a molecular-dynamics trajectory generator. Frame number records position in the requested-amplitude order rather than physical time. The software does not restore side-chain packing, solvent interactions, hydrogen bonding, or full-atom stereochemistry; subsequent rebuilding and refinement remain necessary for full-atom applications.

Article-specific benchmark calculations, manuscript source data, and figure/table reproduction scripts are intentionally excluded from this user release and are distributed in the separate article reproduction package.

## Testing

```bash
python -m pip install ".[test]"
pytest -q
```

The GitHub repository also includes an automated test workflow for Python 3.10–3.13.

## Repository and citation

The source code and fixed software releases are available from:

- Repository: https://github.com/chriszyu9605/BA-ARAP-NMA
- Version 1.0.0 release: https://github.com/chriszyu9605/BA-ARAP-NMA/releases/tag/v1.0.0

If you use BA-ARAP-NMA, please cite the software using the metadata provided in `CITATION.cff`. The citation for the associated article will be added after publication.

## License

MIT License. See `LICENSE`.
