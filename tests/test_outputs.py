import csv
from pathlib import Path
from ba_arap_nma.analysis.geometry import local_geometry_metrics
from ba_arap_nma.core.structure import load_structure
ROOT = Path(__file__).parent

def test_public_metric_set_only():
    structure = load_structure(ROOT / "data/toy_input.pdb")
    metrics = local_geometry_metrics(structure.coords, structure.coords.copy(), structure.segment_ids)
    assert set(metrics) == {
        "mean_squared_adjacent_ca_distance_deviation_A2",
        "adjacent_ca_distance_mae_A",
        "adjacent_distance_outliers_gt_0p5_A",
        "adjacent_distance_outliers_gt_1p0_A",
        "virtual_angle_mae_deg",
        "virtual_angle_outliers_gt_10_deg",
    }
