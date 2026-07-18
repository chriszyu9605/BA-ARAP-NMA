# Changelog

## 1.0.0 — 2026-07-18

- First stable public release.
- Added target-free and two-state command-line workflows.
- Integrated cumulative-overlap, displacement-retention, leading-mode-concentration, and local-distance-change analyses for paired PDB states.
- Standardized PDB, PyMOL, frame-metric, and summary outputs.
- Defined the frame-level squared-distance metric as the parameter-free mean squared adjacent C-alpha distance deviation (A^2), calculated directly from full-precision generated coordinates.
- Documented amplitude-ordered structure sets and clarified that frame indices do not represent physical time.
- Enforced `.pdb` input files at the command-line interface.
