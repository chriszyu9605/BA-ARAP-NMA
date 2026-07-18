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

    security = (root / "SECURITY.md").read_text(encoding="utf-8")
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    license_text = (root / "LICENSE").read_text(encoding="utf-8")

    assert "zhangzhenyu222@mails.ucas.ac.cn" in security
    assert "zhangzhenyu222@mails.ucas.ac.cn" in citation
    assert "zhangzhenyu222@mails.ucas.ac.cn" in pyproject

    for name in ("Zhenyu Zhang", "Dejian Liu", "Haiying Yu", "Luyan Ma"):
        assert name in pyproject
        assert name in license_text

    for family, given in (("Zhang", "Zhenyu"), ("Liu", "Dejian"), ("Yu", "Haiying"), ("Ma", "Luyan")):
        assert f"family-names: {family}" in citation
        assert f"given-names: {given}" in citation
