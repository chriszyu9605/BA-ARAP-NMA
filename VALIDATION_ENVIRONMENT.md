# Validation environment

The release archive was validated with Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, and pytest 9.1.1 on Linux.

The package supports Python 3.10–3.13 through the dependency ranges declared in `pyproject.toml`. `requirements-lock.txt` and `environment.yml` reproduce the exact Python 3.13.5 validation environment; they are not cross-version locks for Python 3.10–3.12.
