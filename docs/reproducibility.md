# Reproducibility

The standard installation uses the dependency ranges declared in `pyproject.toml` and supports Python 3.10–3.13.

`requirements-lock.txt` and `environment.yml` record the exact environment used for final release validation: Python 3.13.5, NumPy 2.3.5, and SciPy 1.17.0. The lock file is not intended as a cross-version lock for Python 3.10–3.12.

The public repository contains the source code, examples, tests, and validation records. Article-specific benchmark calculations, frozen source data, and figure/table reproduction scripts are supplied separately in the article reproduction package.
