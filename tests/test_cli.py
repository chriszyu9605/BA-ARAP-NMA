import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent
CONSOLE = Path(sys.executable).with_name("ba-arap-nma")


def run(*args):
    return subprocess.run([str(CONSOLE), *args], text=True, capture_output=True, check=True)


def test_version():
    assert CONSOLE.exists()
    assert run("version").stdout.strip() == "1.0.0"
    assert run("--version").stdout.strip() == "1.0.0"


def test_generate_cli(tmp_path):
    out = tmp_path / "out"
    result = run(
        "generate", "--pdb", str(ROOT / "data/toy_input.pdb"),
        "--modes", "1", "--orientations", "1",
        "--max-amplitude", "1", "--amplitude-step", "0.5",
        "--output", str(out)
    )
    assert "Completed" in result.stdout
    branch = next(out.glob("ba_arap_nma_mode*"))
    rows = list(csv.DictReader((branch / "frame_metrics.tsv").open(), delimiter="\t"))
    assert len(rows) == 3
    assert (branch / "ca_structures.pdb").exists()
    assert (branch / "view_ca_structures.pml").exists()


def test_two_state_cli(tmp_path):
    out = tmp_path / "two"
    run(
        "two-state", "--input-pdb", str(ROOT / "data/toy_input.pdb"),
        "--target-pdb", str(ROOT / "data/toy_target.pdb"),
        "--input-chain", "A", "--target-chain", "A",
        "--modes", "1", "--orientations", "1",
        "--analysis-modes", "2", "--summary-k", "1",
        "--max-amplitude", "1", "--amplitude-step", "0.5",
        "--min-identity", "0", "--min-coverage", "0",
        "--output", str(out)
    )
    assert (out / "modal_space_analysis.tsv").exists()
    assert (out / "modal_space_summary.tsv").exists()
    assert (out / "best_structures/view_best_structures.pml").exists()
    assert (out / "best_structures/ba_arap_nma_best_ca.pdb").exists()
    assert (out / "best_structures/linear_ca_anm_best_ca.pdb").exists()


def test_fullatom_decorated_output(tmp_path):
    out = tmp_path / "full"
    run(
        "generate", "--pdb", str(ROOT / "data/toy_fullatom_input.pdb"),
        "--chain", "A", "--modes", "1", "--orientations", "1",
        "--max-amplitude", "0.5", "--amplitude-step", "0.5",
        "--output", str(out)
    )
    branch = next(out.glob("ba_arap_nma_mode*"))
    assert (branch / "fullatom_decorated_structures.pdb").exists()
    assert (branch / "view_fullatom_decorated.pml").exists()
    assert (branch / "selected_fullatom_structures/selected_fullatom_structures.tsv").exists()


def test_rejects_non_pdb_input(tmp_path):
    bad = tmp_path / "input.txt"
    bad.write_text((ROOT / "data/toy_input.pdb").read_text())
    result = subprocess.run(
        [str(CONSOLE), "generate", "--pdb", str(bad), "--modes", "1",
         "--orientations", "1", "--max-amplitude", "0",
         "--amplitude-step", "0.5", "--output", str(tmp_path / "bad")],
        text=True, capture_output=True
    )
    assert result.returncode != 0
    assert ".pdb file" in (result.stderr + result.stdout)


def test_output_directory_protection_and_overwrite(tmp_path):
    out = tmp_path / "protected"
    out.mkdir()
    (out / "old.txt").write_text("old")
    base = [
        str(CONSOLE), "generate", "--pdb", str(ROOT / "data/toy_input.pdb"),
        "--modes", "1", "--orientations", "1", "--max-amplitude", "0",
        "--amplitude-step", "0.5", "--output", str(out), "--quiet"
    ]
    rejected = subprocess.run(base, text=True, capture_output=True)
    assert rejected.returncode != 0
    assert "not empty" in (rejected.stderr + rejected.stdout)
    subprocess.run(base + ["--overwrite"], text=True, capture_output=True, check=True)
    assert not (out / "old.txt").exists()
    assert (out / "branch_summary.tsv").exists()
