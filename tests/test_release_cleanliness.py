from pathlib import Path


def test_no_internal_or_unvalidated_content():
    root = Path(__file__).parents[1]
    forbidden = [
        "NOLB" + "_STYLE",
        "NOLB" + "-style",
        "v" + "17",
        "pseudo" + "_dihedral",
        "quality" + "_label",
        "quality" + "_flags",
        "nonlocal" + "_ca",
    ]
    skip = {Path(__file__).resolve()}
    for path in root.rglob("*"):
        if path.resolve() in skip or not path.is_file():
            continue
        if path.suffix.lower() not in {".py", ".md", ".toml", ".txt", ".yml", ".yaml", ".cff", ".sh", ".ps1"}:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        for term in forbidden:
            assert term not in text, (term, path)
