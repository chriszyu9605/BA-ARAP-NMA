# Quick smoke test

This small example confirms that the installed command-line interface can generate structures and complete a two-state analysis in a few seconds.

Linux/macOS:

```bash
bash run_smoke_test.sh
```

Windows PowerShell:

```powershell
./run_smoke_test.ps1
```

The scripts create `output_generate/` and `output_two_state/` in this directory. Existing output directories are replaced explicitly with `--overwrite`.
