# Output files

## Single-state `generate` workflow

Each mode-rank/orientation branch contains:

- `ca_structures.pdb`: amplitude-ordered multi-model C-alpha structures.
- `ca_trace_structures.pdb`: the same structures with sequential `CONECT` records for trace visualization.
- `view_ca_structures.pml`: PyMOL script for viewing all structural states.
- `selected_structures/`: structures at the start, quarter points, midpoint, three-quarter point, and end of the sampled amplitude range.
- `frame_metrics.tsv`: requested amplitude and the local-geometry measurements used in the associated article.
- `branch_summary.json`: branch-level summary.
- `run_config.json`: exact branch settings.
- `fullatom_decorated_structures.pdb`, `selected_fullatom_structures/`, and `view_fullatom_decorated.pml` when the input PDB contains non-C-alpha protein atoms.

The run directory also contains `branch_summary.tsv` and `run_metadata.json`.

## Two-state `two-state` workflow

The two-state workflow creates matched `ba_arap_nma/` and `linear_ca_anm/` branch directories and additionally writes:

- `modal_space_analysis.tsv`: cumulative overlap for each requested K.
- `modal_space_summary.tsv` and `.json`: eta, maximum attainable overlap, rho_K, and s_loc.
- `all_branch_summary.tsv`: all generated branches.
- `best_structures/`: the minimum-RMSD BA-ARAP-NMA and Linear C-alpha ANM structures, the aligned paired second state, a comparison table, and a PyMOL script.

The paired second state is used only for mapping, alignment, modal-space analysis, and retrospective comparison; it does not affect generated BA-ARAP-NMA coordinates.
