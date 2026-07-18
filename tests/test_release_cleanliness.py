from pathlib import Path


def test_public_release_top_level_layout():
    root = Path(__file__).resolve().parents[1]
    required = {
        "README.md", "LICENSE", "CITATION.cff", "pyproject.toml",
        "src", "tests", "docs", "examples", ".github",
    }
    present = {path.name for path in root.iterdir()}
    assert required <= present
