# GitHub release checklist

1. Create a public repository named `BA-ARAP-NMA`. Do not pre-populate it with a README, `.gitignore`, or license because those files are already in this release.
2. Upload the contents of this directory, not the enclosing ZIP file.
3. Confirm that the MIT `LICENSE`, `README.md`, `CITATION.cff`, examples, tests, and `.github/workflows/ci.yml` are visible.
4. Wait for the GitHub Actions test workflow to pass.
5. Create tag `v1.0.0` from this validated initial commit and create a GitHub Release titled `BA-ARAP-NMA 1.0.0`. Attach the validated release ZIP if a separately downloadable archive is desired.
6. After the release URL is known, add the repository URL to `pyproject.toml` under `[project.urls]` and to `CITATION.cff` as `repository-code`, then commit and push this metadata update to `main`. Do not move the `v1.0.0` tag; it should remain attached to the validated release contents.
7. Add the repository or release URL to the manuscript Data Availability Statement.
8. After the associated article is published, add its citation to `CITATION.cff` and the README.
