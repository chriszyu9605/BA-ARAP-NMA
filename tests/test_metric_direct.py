import numpy as np
from ba_arap_nma.analysis.geometry import local_geometry_metrics


def test_mean_squared_adjacent_distance_deviation_is_direct():
    reference = np.array([[0.,0.,0.],[1.,0.,0.],[2.,0.,0.]])
    coords = np.array([[0.,0.,0.],[1.1,0.,0.],[2.3,0.,0.]])
    metrics = local_geometry_metrics(reference, coords, [0,0,0])
    # Deviations are +0.1 and +0.2 A; mean square = (0.01 + 0.04) / 2 = 0.025 A^2.
    assert np.isclose(metrics['mean_squared_adjacent_ca_distance_deviation_A2'], 0.025)
