# Installation

## Standard installation

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
ba-arap-nma --version
```

Windows PowerShell activation:

```powershell
.\.venv\Scripts\Activate.ps1
```

## Exact validated environment

The lock file records the Python 3.13.5 validation environment. Create the environment with Python 3.13.5, then run:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-lock.txt
python -m pip install --no-deps .
python -m pip install "pytest>=8"
pytest -q
```

For Python 3.10–3.12, use the standard installation; `pyproject.toml` selects compatible NumPy and SciPy versions. The validated environment is recorded in `VALIDATION_ENVIRONMENT.md`. PyMOL is optional and is not installed by this package.
