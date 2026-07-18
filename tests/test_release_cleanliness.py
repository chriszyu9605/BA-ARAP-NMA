from pathlib import Path


def test_public_release_top_level_layout():
    root = Path(__file__).resolve().parents[1]
    expected = {
        ".github", ".gitignore", "CHANGELOG.md", "CITATION.cff",
        "CONTRIBUTING.md", "LICENSE", "README.md", "SECURITY.md",
        "THIRD_PARTY_NOTICES.md", "VALIDATION_ENVIRONMENT.md",
        "VALIDATION_REPORT.md", "docs", "environment.yml", "examples",
        "pyproject.toml", "requirements-lock.txt", "requirements.txt",
        "src", "tests",
    }
    ignored = {".git", ".pytest_cache", ".venv", "build", "dist"}
    present = {path.name for path in root.iterdir() if path.name not in ignored}
    assert present == expected
