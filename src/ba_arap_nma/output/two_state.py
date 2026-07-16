from pathlib import Path
from .writers import write_single, write_tsv

METHODS = ("ba_arap_nma", "linear_ca_anm")

def write_best_comparison(outdir, records, target_records, method_best):
    outdir = Path(outdir)
    best_dir = outdir / "best_structures"
    best_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for method in METHODS:
        item = method_best[method]
        filename = f"{method}_best_ca.pdb"
        write_single(best_dir / filename, records, item["coords"])
        row = {key: value for key, value in item.items() if key != "coords"}
        row["file"] = filename
        rows.append(row)
    write_single(
        best_dir / "paired_second_state_ca.pdb",
        target_records,
        method_best["_target"]["coords"],
    )
    write_tsv(best_dir / "best_structure_comparison.tsv", rows)
    (best_dir / "view_best_structures.pml").write_text(
        "reinitialize\n"
        "load paired_second_state_ca.pdb, target\n"
        "load ba_arap_nma_best_ca.pdb, ba\n"
        "load linear_ca_anm_best_ca.pdb, linear\n"
        "hide everything, all\n"
        "set cartoon_trace_atoms, on\n"
        "show cartoon, target or ba or linear\n"
        "show sticks, target or ba or linear\n"
        "color green, target\n"
        "color orange, ba\n"
        "color marine, linear\n"
        "orient target or ba or linear\n"
        "bg_color white\n",
        encoding="utf-8",
    )
