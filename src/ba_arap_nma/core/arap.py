from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import bmat, csc_matrix, coo_matrix, eye, kron
from scipy.sparse.linalg import splu

from . import ca_core as core

Array = NDArray[np.float64]

_ARAP_KKT_CACHE: dict[tuple, tuple[object, csc_matrix, list[tuple[int,int,float]], list[list[int]]]] = {}


@dataclass
class ProjectionAudit:
    tangent_cosine: float
    correction_rms_A: float
    correction_to_predictor_ratio: float
    arap_iterations: int
    component_count: int
    component_gauge: str
    rotation_policy: str


def connected_components_from_edges(n: int, edges: Sequence[tuple[int, int, float]]) -> list[list[int]]:
    parent = list(range(n))
    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    for i, j, _ in edges:
        union(i, j)
    comps: dict[int, list[int]] = {}
    for i in range(n):
        comps.setdefault(find(i), []).append(i)
    return list(comps.values())


def _centroid_constraint_matrix(n: int, comps: list[list[int]], target_coords: Array) -> tuple[csc_matrix, Array]:
    rows = []
    cols = []
    vals = []
    rhs = []
    row = 0
    for comp in comps:
        if not comp:
            continue
        w = 1.0 / len(comp)
        centroid = target_coords[comp].mean(axis=0)
        for axis in range(3):
            for idx in comp:
                rows.append(row)
                cols.append(3 * idx + axis)
                vals.append(w)
            rhs.append(float(centroid[axis]))
            row += 1
    mat = coo_matrix((vals, (rows, cols)), shape=(row, 3 * n)).tocsc()
    return mat, np.asarray(rhs, dtype=float)


def arap_tangent_preserving_project(
    reference: Array,
    predictor: Array,
    segment_ids: NDArray[np.int64],
    graph: str = "seq2",
    component_gauge_target: str = "predictor",
) -> tuple[Array, ProjectionAudit]:
    """Component-gauged tangent-preserving ARAP projection.

    Apply a start-referenced co-rotational ARAP correction.

    Independent centroid gauges remove arbitrary translation of disconnected
    local C-alpha segments without introducing interactions between them.
    """
    n = len(reference)
    if predictor.shape != reference.shape:
        raise ValueError("predictor/reference shape mismatch")
    edges = core.build_local_edges(segment_ids, graph=graph)
    if not edges:
        raise ValueError("ARAP local graph has no edges")

    neighbors = core._neighbors_from_edges(n, edges)
    rotations = core._best_rotations(reference, predictor, neighbors)
    identity = np.eye(3)
    remainders: list[Array] = []
    for r in rotations:
        omega = core._rotation_log_skew(r)
        remainders.append(r - identity - omega)

    edge_targets: list[Array] = []
    for i, j, _w in edges:
        ref_edge = reference[i] - reference[j]
        pred_edge = predictor[i] - predictor[j]
        nonlinear_rotation = 0.5 * (remainders[i] + remainders[j]) @ ref_edge
        edge_targets.append(pred_edge + nonlinear_rotation)

    # The left-hand KKT matrix depends only on the local graph and component gauge
    # topology, not on the frame coordinates.  Cache the sparse LU factorization
    # per protein/segment pattern to avoid refactorizing it for every amplitude frame.
    cache_key = (n, tuple(int(x) for x in segment_ids.tolist()), graph)
    cached = _ARAP_KKT_CACHE.get(cache_key)
    if cached is None:
        lap = core._laplacian_from_edges(n, edges)
        l3 = kron(lap, eye(3, format="csr"), format="csr")
        comps = connected_components_from_edges(n, edges)
        # Matrix pattern only; RHS uses the current gauge target below.
        a_sparse, _rhs_dummy = _centroid_constraint_matrix(n, comps, reference)
        zero = csc_matrix((a_sparse.shape[0], a_sparse.shape[0]))
        kkt = bmat([[l3.tocsc(), a_sparse.T], [a_sparse, zero]], format="csc")
        kkt = kkt + eye(kkt.shape[0], format="csc") * 1e-12
        lu = splu(kkt)
        cached = (lu, a_sparse, edges, comps)
        _ARAP_KKT_CACHE[cache_key] = cached
    lu, a_sparse, _edges_cached, comps = cached
    if component_gauge_target == "reference":
        gauge_target = reference
    elif component_gauge_target == "predictor":
        gauge_target = predictor
    else:
        raise ValueError("component_gauge_target must be 'reference' or 'predictor'")

    _a_unused, rhs_gauge = _centroid_constraint_matrix(n, comps, gauge_target)
    rhs = np.concatenate([core._rhs_from_edge_targets(n, edges, edge_targets), rhs_gauge])
    solution = lu.solve(rhs)
    z = solution[: 3 * n].reshape(n, 3)

    # Remove residual global rigid drift relative to the previous conformation.
    z, _, _ = core.kabsch_align(z, reference)

    raw = predictor.reshape(-1) - reference.reshape(-1)
    actual = z.reshape(-1) - reference.reshape(-1)
    raw_norm = float(np.linalg.norm(raw))
    actual_norm = float(np.linalg.norm(actual))
    tangent_cos = float((raw @ actual) / max(raw_norm * actual_norm, 1e-12))
    correction = z - predictor
    correction_rms = float(np.sqrt(np.mean(np.sum(correction * correction, axis=1))))
    predictor_rms = core.ca_rmsd(predictor, reference, align=False)
    return z, ProjectionAudit(
        tangent_cosine=tangent_cos,
        correction_rms_A=correction_rms,
        correction_to_predictor_ratio=float(correction_rms / max(predictor_rms, 1e-12)),
        arap_iterations=1,
        component_count=len(comps),
        component_gauge=f"centroid_per_component_to_{component_gauge_target}",
        rotation_policy="tangent_preserving_second_order_component_gauged",
    )
