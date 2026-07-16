from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence
import numpy as np
from numpy.typing import NDArray
from scipy.sparse import coo_matrix, csr_matrix, eye
from scipy.sparse.linalg import eigsh
from scipy.linalg import qr as dense_qr, eigh as dense_eigh
from scipy.spatial import cKDTree

Array = NDArray[np.float64]

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "C", "PYL": "K",
}
@dataclass(frozen=True)
class CARecord:
    chain: str
    resseq: int
    icode: str
    resname: str
    coord: Array
    serial: int
    source_index: int

    @property
    def aa(self) -> str:
        return AA3_TO_1.get(self.resname.upper(), "X")

    @property
    def label(self) -> str:
        return f"{self.chain}:{self.resname}{self.resseq}{self.icode.strip()}"

@dataclass(frozen=True)
class ParseAudit:
    accepted_atom_ca_count: int
    accepted_modified_residue_ca_count: int
    excluded_hetero_ca_count: int
    excluded_unknown_residue_ca_count: int

@dataclass
class PairData:
    start_records: list[CARecord]
    target_records: list[CARecord]
    start_coords: Array
    target_coords_aligned: Array
    sequence_identity: float
    aligned_count: int
    start_chain_total: int
    target_chain_total: int
    segment_ids: NDArray[np.int64]
    start_parse_audit: ParseAudit
    target_parse_audit: ParseAudit

@dataclass
class Modes:
    eigenvalues: Array
    vectors: Array  # shape (3N, M), Euclidean-orthonormal
    zero_mode_count: int
    anm_edge_count: int
    connected_components: int
    rigid_basis_rank: int
    rigid_projection_max_abs: float
    constraint_scheme: str = "NONE"
    constraint_count: int = 0
    constraint_rank: int = 0
    constrained_nullity: int = 0
    constraint_projection_max_abs: float = 0.0

def parse_ca_records_with_audit(path: str | Path, chain: str | None = None) -> tuple[list[CARecord], ParseAudit]:
    """Read protein C-alpha atoms while excluding ligands and metal ions.

    Standard protein residues are accepted from ATOM records.  A small set of
    genetically encoded/commonly substituted amino acids (MSE, SEC, PYL) is
    also accepted when represented as HETATM.  All other HETATM records named
    ``CA`` (for example calcium ions or ligand atoms) are explicitly excluded.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    out: list[CARecord] = []
    seen: set[tuple[str, int, str]] = set()
    source_index = 0
    accepted_atom = 0
    accepted_modified = 0
    excluded_hetero = 0
    excluded_unknown = 0
    allowed_modified = {"MSE", "SEC", "PYL"}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record_type = line[:6]
            if record_type not in {"ATOM  ", "HETATM"}:
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            altloc = line[16:17]
            if altloc not in (" ", "A", "1"):
                continue
            this_chain = line[21:22].strip() or "_"
            if chain is not None and this_chain != chain:
                continue
            resname = line[17:20].strip().upper()
            if record_type == "HETATM" and resname not in allowed_modified:
                excluded_hetero += 1
                continue
            if resname not in AA3_TO_1:
                excluded_unknown += 1
                continue
            try:
                resseq = int(line[22:26])
                icode = line[26:27]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                serial = int(line[6:11])
            except ValueError:
                continue
            key = (this_chain, resseq, icode)
            if key in seen:
                continue
            seen.add(key)
            out.append(CARecord(
                chain=this_chain,
                resseq=resseq,
                icode=icode,
                resname=resname,
                coord=np.array([x, y, z], dtype=float),
                serial=serial,
                source_index=source_index,
            ))
            source_index += 1
            if record_type == "ATOM  ":
                accepted_atom += 1
            else:
                accepted_modified += 1
    if not out:
        raise ValueError(f"No protein C-alpha atoms found in {path} for chain={chain!r}")
    return out, ParseAudit(
        accepted_atom_ca_count=accepted_atom,
        accepted_modified_residue_ca_count=accepted_modified,
        excluded_hetero_ca_count=excluded_hetero,
        excluded_unknown_residue_ca_count=excluded_unknown,
    )

def parse_ca_records(path: str | Path, chain: str | None = None) -> list[CARecord]:
    records, _ = parse_ca_records_with_audit(path, chain)
    return records

def needleman_wunsch(seq_a: str, seq_b: str, match: int = 2, mismatch: int = -1, gap: int = -2) -> list[tuple[int | None, int | None]]:
    n, m = len(seq_a), len(seq_b)
    score = np.zeros((n + 1, m + 1), dtype=np.int32)
    trace = np.zeros((n + 1, m + 1), dtype=np.int8)  # 0 diag, 1 up, 2 left
    score[:, 0] = np.arange(n + 1) * gap
    score[0, :] = np.arange(m + 1) * gap
    trace[1:, 0] = 1
    trace[0, 1:] = 2
    for i in range(1, n + 1):
        ai = seq_a[i - 1]
        for j in range(1, m + 1):
            diag = score[i - 1, j - 1] + (match if ai == seq_b[j - 1] else mismatch)
            up = score[i - 1, j] + gap
            left = score[i, j - 1] + gap
            best = max(diag, up, left)
            score[i, j] = best
            trace[i, j] = 0 if best == diag else (1 if best == up else 2)
    pairs: list[tuple[int | None, int | None]] = []
    i, j = n, m
    while i > 0 or j > 0:
        direction = trace[i, j]
        if i > 0 and j > 0 and direction == 0:
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or direction == 1):
            pairs.append((i - 1, None))
            i -= 1
        else:
            pairs.append((None, j - 1))
            j -= 1
    pairs.reverse()
    return pairs

def kabsch_align(moving: Array, reference: Array) -> tuple[Array, Array, Array]:
    if moving.shape != reference.shape:
        raise ValueError("moving/reference shape mismatch")
    cm = moving.mean(axis=0)
    cr = reference.mean(axis=0)
    x = moving - cm
    y = reference - cr
    cov = x.T @ y
    u, _, vt = np.linalg.svd(cov)
    rot = vt.T @ u.T
    if np.linalg.det(rot) < 0:
        vt[-1, :] *= -1
        rot = vt.T @ u.T
    trans = cr - cm @ rot.T
    aligned = moving @ rot.T + trans
    return aligned, rot, trans

def prepare_pair(start_pdb: str | Path, target_pdb: str | Path, start_chain: str, target_chain: str, min_identity: float = 0.85, min_coverage: float = 0.80) -> PairData:
    s_all, s_audit = parse_ca_records_with_audit(start_pdb, start_chain)
    t_all, t_audit = parse_ca_records_with_audit(target_pdb, target_chain)
    seq_s = "".join(r.aa for r in s_all)
    seq_t = "".join(r.aa for r in t_all)
    alignment = needleman_wunsch(seq_s, seq_t)
    mapped: list[tuple[int, int]] = [(i, j) for i, j in alignment if i is not None and j is not None]
    if not mapped:
        raise ValueError("No aligned residues")
    identical = sum(seq_s[i] == seq_t[j] for i, j in mapped)
    identity = identical / len(mapped)
    coverage = len(mapped) / min(len(s_all), len(t_all))
    if identity < min_identity:
        raise ValueError(f"Sequence identity too low: {identity:.3f} < {min_identity}")
    if coverage < min_coverage:
        raise ValueError(f"Alignment coverage too low: {coverage:.3f} < {min_coverage}")
    s_records = [s_all[i] for i, _ in mapped]
    t_records = [t_all[j] for _, j in mapped]
    s = np.stack([r.coord for r in s_records])
    t = np.stack([r.coord for r in t_records])
    t_aligned, _, _ = kabsch_align(t, s)

    # Break local backbone segments across missing residue-number gaps.
    segment_ids = np.zeros(len(s_records), dtype=np.int64)
    seg = 0
    for k in range(1, len(s_records)):
        prev, cur = s_records[k - 1], s_records[k]
        if cur.chain != prev.chain or (cur.resseq - prev.resseq) > 1 or (cur.resseq - prev.resseq) < 0:
            seg += 1
        segment_ids[k] = seg
    return PairData(
        start_records=s_records,
        target_records=t_records,
        start_coords=s,
        target_coords_aligned=t_aligned,
        sequence_identity=identity,
        aligned_count=len(mapped),
        start_chain_total=len(s_all),
        target_chain_total=len(t_all),
        segment_ids=segment_ids,
        start_parse_audit=s_audit,
        target_parse_audit=t_audit,
    )

def _count_components(n: int, edges: Sequence[tuple[int, int]]) -> int:
    parent = list(range(n))
    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    for i, j in edges:
        union(i, j)
    return len({find(i) for i in range(n)})

def build_anm_hessian(coords: Array, cutoff_A: float = 15.0, gamma: float = 1.0) -> tuple[csr_matrix, int, int]:
    n = len(coords)
    tree = cKDTree(coords)
    pairs = sorted(tree.query_pairs(cutoff_A))
    if not pairs:
        raise ValueError("ANM graph has no edges")
    components = _count_components(n, pairs)
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    diag = np.zeros((n, 3, 3), dtype=float)
    for i, j in pairs:
        vec = coords[j] - coords[i]
        r2 = float(vec @ vec)
        if r2 < 1e-12:
            continue
        block = gamma * np.outer(vec, vec) / r2
        diag[i] += block
        diag[j] += block
        for a in range(3):
            for b in range(3):
                value = -block[a, b]
                if value != 0.0:
                    rows.extend([3 * i + a, 3 * j + a])
                    cols.extend([3 * j + b, 3 * i + b])
                    data.extend([value, value])
    for i in range(n):
        for a in range(3):
            for b in range(3):
                value = diag[i, a, b]
                if value != 0.0:
                    rows.append(3 * i + a)
                    cols.append(3 * i + b)
                    data.append(value)
    h = coo_matrix((data, (rows, cols)), shape=(3 * n, 3 * n)).tocsr()
    h = (h + h.T) * 0.5
    return h, len(pairs), components

def rigid_body_basis(coords: Array) -> Array:
    """Return an orthonormal basis for global translations and rotations.

    For a non-degenerate 3D protein this basis has rank six.  The basis is
    constructed analytically rather than inferred from small eigenvalues.
    """
    n = len(coords)
    centered = coords - coords.mean(axis=0)
    cols: list[Array] = []
    for axis in range(3):
        col = np.zeros((n, 3), dtype=float)
        col[:, axis] = 1.0
        cols.append(col.reshape(-1))
    axes = np.eye(3)
    for axis in axes:
        # Infinitesimal rigid rotation: omega x r.
        cols.append(np.cross(np.broadcast_to(axis, centered.shape), centered).reshape(-1))
    raw = np.column_stack(cols)
    u, singular, _ = np.linalg.svd(raw, full_matrices=False)
    tol = max(raw.shape) * np.finfo(float).eps * max(float(singular[0]), 1.0)
    rank = int(np.sum(singular > tol))
    return np.asarray(u[:, :rank], dtype=float)

def compute_low_modes(coords: Array, n_modes: int = 10, cutoff_A: float = 15.0, eig_tol: float = 1e-7) -> Modes:
    """Compute low ANM modes with explicit removal of global rigid motions.

    The sparse eigensolver first obtains a small low-eigenvalue subspace.
    The six analytical rigid-body vectors are then projected out, followed
    by a Rayleigh--Ritz diagonalization in the projected subspace.  This is
    more reliable than deciding whether five or six modes are zero from an
    eigenvalue threshold.
    """
    h, edge_count, components = build_anm_hessian(coords, cutoff_A=cutoff_A)
    if components != 1:
        raise ValueError(f"ANM graph is disconnected: {components} components")
    dof = h.shape[0]
    rigid = rigid_body_basis(coords)
    rigid_rank = int(rigid.shape[1])
    if rigid_rank != 6:
        raise ValueError(f"Expected six rigid-body vectors, obtained rank {rigid_rank}")

    raw_count = min(dof, n_modes + rigid_rank + 12)
    if dof <= raw_count + 1:
        vals, vecs = np.linalg.eigh(h.toarray())
        order = np.argsort(vals)
        vals = np.asarray(vals[order], dtype=float)
        vecs = np.asarray(vecs[:, order], dtype=float)
        vecs = vecs[:, :raw_count]
    else:
        shifted = h + eye(dof, format="csr") * 1e-10
        vals, vecs = eigsh(
            shifted,
            k=min(dof - 2, raw_count),
            which="SM",
            tol=eig_tol,
            maxiter=max(20000, dof * 30),
        )
        vals = np.asarray(vals - 1e-10, dtype=float)
        order = np.argsort(vals)
        vals = vals[order]
        vecs = np.asarray(vecs[:, order], dtype=float)

    projected = vecs - rigid @ (rigid.T @ vecs)
    u, singular, _ = np.linalg.svd(projected, full_matrices=False)
    sv_tol = max(projected.shape) * np.finfo(float).eps * max(float(singular[0]), 1.0)
    keep = singular > sv_tol
    basis = np.asarray(u[:, keep], dtype=float)
    if basis.shape[1] < n_modes:
        raise ValueError(
            f"Only {basis.shape[1]} non-rigid vectors available; requested {n_modes}. "
            "Increase the raw eigenspace size."
        )

    reduced = basis.T @ (h @ basis)
    reduced = 0.5 * (reduced + reduced.T)
    rvals, rvecs = np.linalg.eigh(reduced)
    order = np.argsort(rvals)
    rvals = np.asarray(rvals[order], dtype=float)
    rvecs = np.asarray(rvecs[:, order], dtype=float)
    threshold = max(1e-10, float(np.max(np.abs(rvals))) * 1e-10)
    positive = np.where(rvals > threshold)[0]
    if len(positive) < n_modes:
        raise ValueError(f"Only {len(positive)} positive projected modes found; requested {n_modes}")
    idx = positive[:n_modes]
    selected_vecs = basis @ rvecs[:, idx]
    # Re-project and orthonormalize once more to suppress round-off.
    selected_vecs = selected_vecs - rigid @ (rigid.T @ selected_vecs)
    selected_vecs, _ = np.linalg.qr(selected_vecs)
    selected_vecs = selected_vecs[:, :n_modes]
    selected_vals = np.array([float(v @ (h @ v)) for v in selected_vecs.T], dtype=float)
    eig_order = np.argsort(selected_vals)
    selected_vals = selected_vals[eig_order]
    selected_vecs = selected_vecs[:, eig_order]
    rigid_residual = float(np.max(np.abs(rigid.T @ selected_vecs)))
    return Modes(
        eigenvalues=selected_vals,
        vectors=selected_vecs,
        zero_mode_count=rigid_rank,
        anm_edge_count=edge_count,
        connected_components=components,
        rigid_basis_rank=rigid_rank,
        rigid_projection_max_abs=rigid_residual,
        constraint_scheme="NONE",
        constraint_count=0,
        constraint_rank=0,
        constrained_nullity=3 * len(coords) - rigid_rank,
        constraint_projection_max_abs=0.0,
    )

def build_local_distance_constraint_jacobian(
    coords: Array,
    segment_ids: NDArray[np.int64],
    scheme: str,
) -> tuple[Array, list[tuple[int, int]]]:
    """Build exact first-order C-alpha distance constraints.

    C1 constrains only adjacent C-alpha distances (i,i+1).  C2 adds the
    (i,i+2) distances, which also constrains the virtual backbone angle at
    first order.  Rows use unit edge directions, so all constraints share the
    same physical scale and no empirical edge weights are introduced.
    """
    scheme_u = scheme.upper()
    if scheme_u not in {"C1", "C2"}:
        raise ValueError(f"Unknown local constraint scheme {scheme!r}; expected C1 or C2")
    max_gap = 1 if scheme_u == "C1" else 2
    n = len(coords)
    rows: list[Array] = []
    pairs: list[tuple[int, int]] = []
    for i in range(n):
        for gap in range(1, max_gap + 1):
            j = i + gap
            if j >= n:
                break
            if segment_ids[i] != segment_ids[j]:
                break
            edge = coords[i] - coords[j]
            length = float(np.linalg.norm(edge))
            if length < 1e-10:
                raise ValueError(f"Degenerate C-alpha constraint edge ({i},{j})")
            unit = edge / length
            row = np.zeros(3 * n, dtype=float)
            row[3 * i:3 * i + 3] = unit
            row[3 * j:3 * j + 3] = -unit
            rows.append(row)
            pairs.append((i, j))
    if not rows:
        raise ValueError("No local distance constraints could be built")
    return np.vstack(rows), pairs

def _orthogonal_complement_basis(
    row_constraints: Array,
    dof: int,
    rank_rtol: float = 1e-10,
) -> tuple[Array, int]:
    """Return an orthonormal basis of the null space of row_constraints.

    A rank-revealing QR decomposition is used instead of a penalty term.  The
    resulting basis is exact up to numerical precision and introduces no
    stiffness weight or protein-specific parameter.
    """
    if row_constraints.ndim != 2 or row_constraints.shape[1] != dof:
        raise ValueError("constraint matrix shape mismatch")
    # QR of A^T: the leading Q columns span range(A^T), the remainder spans
    # null(A).  Column pivoting makes the numerical rank test robust to any
    # redundant local constraints.
    q, r, _piv = dense_qr(row_constraints.T, mode="full", pivoting=True)
    diag = np.abs(np.diag(r))
    scale = max(float(diag.max()) if diag.size else 0.0, 1.0)
    rank = int(np.sum(diag > rank_rtol * scale))
    null_basis = np.asarray(q[:, rank:], dtype=float)
    if null_basis.shape[1] == 0:
        raise ValueError("Local constraints leave no conformational degrees of freedom")
    return null_basis, rank

def compute_constraint_compatible_modes(
    coords: Array,
    segment_ids: NDArray[np.int64],
    scheme: str,
    n_modes: int = 10,
    cutoff_A: float = 15.0,
    rank_rtol: float = 1e-10,
) -> Modes:
    """Compute ANM modes inside an exact local-geometry tangent space.

    The standard C-alpha ANM Hessian supplies the elastic ordering.  Local
    distance Jacobians and the six analytical rigid-body vectors are imposed
    as exact linear constraints.  The Hessian is then diagonalized in their
    common null space.  No penalty stiffness or ARAP/NMA mixing coefficient is
    used.
    """
    h, edge_count, components = build_anm_hessian(coords, cutoff_A=cutoff_A)
    if components != 1:
        raise ValueError(f"ANM graph is disconnected: {components} components")
    dof = h.shape[0]
    rigid = rigid_body_basis(coords)
    if rigid.shape[1] != 6:
        raise ValueError(f"Expected six rigid-body vectors, obtained rank {rigid.shape[1]}")
    j_local, _pairs = build_local_distance_constraint_jacobian(coords, segment_ids, scheme)
    combined = np.vstack([j_local, rigid.T])
    null_basis, combined_rank = _orthogonal_complement_basis(combined, dof, rank_rtol=rank_rtol)
    if null_basis.shape[1] < n_modes:
        raise ValueError(
            f"Constraint scheme {scheme} leaves only {null_basis.shape[1]} non-rigid DOF; "
            f"requested {n_modes} modes"
        )

    reduced = null_basis.T @ (h @ null_basis)
    reduced = 0.5 * (reduced + reduced.T)
    # Only the lowest portion is needed.  The constraint null space has no
    # analytical rigid modes because those were included in `combined`.
    upper = min(reduced.shape[0] - 1, n_modes + 8)
    vals, vecs = dense_eigh(reduced, subset_by_index=[0, upper], driver="evr")
    vals = np.asarray(vals, dtype=float)
    vecs = np.asarray(vecs, dtype=float)
    threshold = max(1e-10, float(np.max(np.abs(vals))) * 1e-10)
    positive = np.where(vals > threshold)[0]
    if len(positive) < n_modes:
        # Fall back to the full reduced eigendecomposition only when the
        # requested window contains unexpected near-zero numerical modes.
        vals, vecs = np.linalg.eigh(reduced)
        vals = np.asarray(vals, dtype=float)
        vecs = np.asarray(vecs, dtype=float)
        positive = np.where(vals > max(1e-10, float(np.max(np.abs(vals))) * 1e-10))[0]
    if len(positive) < n_modes:
        raise ValueError(f"Only {len(positive)} positive constrained modes found; requested {n_modes}")
    idx = positive[:n_modes]
    selected_vals = vals[idx]
    selected_vecs = null_basis @ vecs[:, idx]
    # Numerical cleanup without changing the subspace.
    selected_vecs, _ = np.linalg.qr(selected_vecs)
    selected_vecs = selected_vecs[:, :n_modes]
    selected_vals = np.array([float(v @ (h @ v)) for v in selected_vecs.T], dtype=float)
    order = np.argsort(selected_vals)
    selected_vals = selected_vals[order]
    selected_vecs = selected_vecs[:, order]
    rigid_residual = float(np.max(np.abs(rigid.T @ selected_vecs)))
    local_residual = float(np.max(np.abs(j_local @ selected_vecs)))
    return Modes(
        eigenvalues=selected_vals,
        vectors=selected_vecs,
        zero_mode_count=6,
        anm_edge_count=edge_count,
        connected_components=components,
        rigid_basis_rank=6,
        rigid_projection_max_abs=rigid_residual,
        constraint_scheme=scheme.upper(),
        constraint_count=int(j_local.shape[0]),
        constraint_rank=int(combined_rank - 6),
        constrained_nullity=int(null_basis.shape[1]),
        constraint_projection_max_abs=local_residual,
    )

def build_local_edges(segment_ids: NDArray[np.int64], graph: str = "seq2") -> list[tuple[int, int, float]]:
    if graph not in {"seq1", "seq2", "seq4"}:
        raise ValueError(f"Unknown ARAP graph {graph}")
    max_gap = {"seq1": 1, "seq2": 2, "seq4": 4}[graph]
    n = len(segment_ids)
    edges: list[tuple[int, int, float]] = []
    for i in range(n):
        for gap in range(1, max_gap + 1):
            j = i + gap
            if j >= n:
                break
            if segment_ids[i] != segment_ids[j]:
                break
            edges.append((i, j, 1.0))
    return edges

def _neighbors_from_edges(n: int, edges: Sequence[tuple[int, int, float]]) -> list[list[tuple[int, float]]]:
    neigh: list[list[tuple[int, float]]] = [[] for _ in range(n)]
    for i, j, w in edges:
        neigh[i].append((j, w))
        neigh[j].append((i, w))
    return neigh

def _best_rotations(reference: Array, deformed: Array, neighbors: list[list[tuple[int, float]]]) -> list[Array]:
    rotations: list[Array] = []
    for i, nbrs in enumerate(neighbors):
        if len(nbrs) < 2:
            rotations.append(np.eye(3))
            continue
        cov = np.zeros((3, 3), dtype=float)
        for j, w in nbrs:
            p = reference[j] - reference[i]
            q = deformed[j] - deformed[i]
            cov += w * np.outer(p, q)
        u, _, vt = np.linalg.svd(cov)
        r = vt.T @ u.T
        if np.linalg.det(r) < 0:
            vt[-1, :] *= -1
            r = vt.T @ u.T
        rotations.append(r)
    return rotations

def _laplacian_from_edges(n: int, edges: Sequence[tuple[int, int, float]]) -> csr_matrix:
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    diag = np.zeros(n, dtype=float)
    for i, j, w in edges:
        diag[i] += w
        diag[j] += w
        rows.extend([i, j])
        cols.extend([j, i])
        vals.extend([-w, -w])
    for i, d in enumerate(diag):
        rows.append(i)
        cols.append(i)
        vals.append(float(d))
    return coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()

def _rotation_log_skew(rotation: Array) -> Array:
    """Return the principal matrix logarithm of a proper 3-D rotation.

    The result is a 3x3 skew-symmetric matrix Omega such that exp(Omega)=R
    for rotation angles in [0, pi].  A series-safe expression is used near
    zero.  This routine is intentionally local and does not unwrap rotations
    across large steps; the v04 method is an incremental continuation.
    """
    r = np.asarray(rotation, dtype=float)
    cos_theta = float(np.clip((np.trace(r) - 1.0) * 0.5, -1.0, 1.0))
    theta = float(np.arccos(cos_theta))
    skew_part = 0.5 * (r - r.T)
    if theta < 1e-6:
        return skew_part
    sin_theta = float(np.sin(theta))
    if abs(sin_theta) < 1e-7:
        # Near pi the principal axis is numerically delicate.  The predictor
        # steps used by this method are small, so reaching this branch is a
        # diagnostic failure rather than a place for a hidden heuristic.
        raise ValueError(f"Local predictor rotation too large for incremental log map: {theta:.6g} rad")
    return (theta / sin_theta) * skew_part

def _rhs_from_edge_targets(
    n: int,
    edges: Sequence[tuple[int, int, float]],
    targets: Sequence[Array],
) -> Array:
    if len(edges) != len(targets):
        raise ValueError("edge/target length mismatch")
    b = np.zeros((n, 3), dtype=float)
    for (i, j, w), target in zip(edges, targets):
        term = float(w) * np.asarray(target, dtype=float)
        b[i] += term
        b[j] -= term
    return b.reshape(-1)

def ca_rmsd(a: Array, b: Array, align: bool = True) -> float:
    if align:
        a, _, _ = kabsch_align(a, b)
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))

def transition_coverage(initial: Array, predicted: Array, target: Array) -> float:
    initial_rmsd = ca_rmsd(initial, target, align=True)
    final_rmsd = ca_rmsd(predicted, target, align=True)
    if initial_rmsd < 1e-12:
        return 1.0
    return float((initial_rmsd - final_rmsd) / initial_rmsd)
