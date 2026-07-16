# BA-ARAP-NMA 1.0.0 validation report

The release was validated from a clean Python virtual environment.

- Source installation with `python -m pip install .`: passed.
- Installed command `ba-arap-nma`: passed.
- `ba-arap-nma --version` and `ba-arap-nma version`: passed.
- Target-free `generate` workflow: passed.
- Two-state workflow and modal-space output: passed.
- PDB-only input validation: passed.
- Non-empty output-directory protection: passed.
- Explicit `--overwrite` replacement: passed.
- Progress output and `--quiet`: passed.
- C2 constraint residual test: passed.
- Six-rigid-body-direction removal test: passed.
- 81-frame amplitude-grid test: passed.
- Target-independence test: passed.
- Decorated full-atom PDB/PyMOL output test: passed.
- Public metric set and output-name checks: passed.
- Quick smoke-test example: passed.
- HCV helicase example: archived expected outputs and an article-preparation validation record are supplied. This real-protein example is computationally intensive and runtime depends strongly on the CPU and linear-algebra backend; it is not part of the automated CI test suite.
- Automated test suite: 10 passed.

Validated environment: Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, Linux.
