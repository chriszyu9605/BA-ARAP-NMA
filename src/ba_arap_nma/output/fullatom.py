from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence
import numpy as np
from numpy.typing import NDArray

from ba_arap_nma.core.ca_core import CARecord, kabsch_align

Array = NDArray[np.float64]

@dataclass
class AtomRecord:
    line: str
    record_type: str
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    coord: Array
    element: str

    @property
    def key(self) -> tuple[str,int,str]:
        return (self.chain, self.resseq, self.icode)


def parse_atom_records(path: str | Path, chain: str | None = None, keep_hetero: bool = False) -> list[AtomRecord]:
    path = Path(path)
    atoms: list[AtomRecord] = []
    with path.open('r', encoding='utf-8', errors='replace') as f:
        for line in f:
            rec = line[:6]
            if rec not in {'ATOM  ', 'HETATM'}:
                continue
            if rec == 'HETATM' and not keep_hetero:
                continue
            altloc = line[16:17]
            if altloc not in (' ', 'A', '1'):
                continue
            ch = line[21:22].strip() or '_'
            if chain is not None and ch != chain:
                continue
            try:
                serial = int(line[6:11])
                name = line[12:16]
                resname = line[17:20].strip().upper()
                resseq = int(line[22:26])
                icode = line[26:27]
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except Exception:
                continue
            element = line[76:78].strip() if len(line) >= 78 else ''
            atoms.append(AtomRecord(line=line.rstrip('\n'), record_type=rec.strip(), serial=serial, name=name, altloc=altloc, resname=resname, chain=ch, resseq=resseq, icode=icode, coord=np.array([x,y,z], dtype=float), element=element))
    return atoms


def _format_atom_line(serial: int, atom: AtomRecord, xyz: Array) -> str:
    # Preserve atom/residue identity, occupancy/B-factor if available, update coordinates.
    line = atom.line
    if len(line) < 80:
        line = line.ljust(80)
    rec = (atom.record_type if atom.record_type in {'ATOM','HETATM'} else 'ATOM').ljust(6)
    prefix = f"{rec}{serial:5d} {atom.name:4s}{atom.altloc:1s}{atom.resname:>3s} {(' ' if atom.chain=='_' else atom.chain):1s}{atom.resseq:4d}{atom.icode:1s}   "
    suffix = line[54:]
    return f"{prefix}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{suffix}\n"


def _ca_key_to_index(records: Sequence[CARecord]) -> dict[tuple[str,int,str], int]:
    return {(r.chain, r.resseq, r.icode): i for i, r in enumerate(records)}


def _segment_id_by_records(records: Sequence[CARecord]) -> list[int]:
    seg=[0]*len(records)
    cur=0
    for i in range(1,len(records)):
        a,b=records[i-1],records[i]
        if b.chain!=a.chain or (b.resseq-a.resseq)>1 or (b.resseq-a.resseq)<0:
            cur+=1
        seg[i]=cur
    return seg


def _local_indices(i: int, records: Sequence[CARecord], segs: Sequence[int], radius: int = 2) -> list[int]:
    lo=max(0,i-radius); hi=min(len(records), i+radius+1)
    idx=[j for j in range(lo,hi) if segs[j]==segs[i]]
    if len(idx) < 3:
        idx=[j for j in range(max(0,i-1), min(len(records),i+2)) if segs[j]==segs[i]]
    return idx


def decorate_frame_atoms(ca_records: Sequence[CARecord], start_ca: Array, frame_ca: Array, atom_records: Sequence[AtomRecord], keep_unmapped: bool = False) -> list[tuple[AtomRecord, Array]]:
    """Return updated atom coordinates by residue-local Kabsch decoration.

    The computed model remains Cα-driven.  Full-atom coordinates are generated
    by carrying each residue's original atoms with the local Cα frame.  Side-chain
    rotamers are not resampled; the output is a user-facing structural hypothesis,
    not an atomistically minimized model.
    """
    key_to_idx=_ca_key_to_index(ca_records)
    segs=_segment_id_by_records(ca_records)
    rot_cache: dict[int, tuple[Array, Array, Array, Array]] = {}
    out: list[tuple[AtomRecord, Array]] = []
    for atom in atom_records:
        idx = key_to_idx.get(atom.key)
        if idx is None:
            if keep_unmapped:
                out.append((atom, atom.coord.copy()))
            continue
        if idx not in rot_cache:
            inds=_local_indices(idx, ca_records, segs, radius=2)
            if len(inds) >= 3:
                ref=start_ca[inds]
                mov=frame_ca[inds]
                # kabsch_align returns aligned=moving@rot.T+trans to reference; we need start->frame.
                # So align start-local to frame-local.
                _, rot, trans = kabsch_align(ref, mov)
            elif len(inds) >= 1:
                rot=np.eye(3)
                trans=frame_ca[idx]-start_ca[idx]
            else:
                rot=np.eye(3); trans=np.zeros(3)
            rot_cache[idx]=(rot,trans,start_ca[idx],frame_ca[idx])
        rot, trans, ca0, cat = rot_cache[idx]
        # Use local rigid transform inferred by Kabsch.  Applying as p@rot.T+trans.
        xyz = atom.coord @ rot.T + trans
        # Stabilize CA atom exactly at generated Cα coordinate.
        if atom.name.strip() == 'CA':
            xyz = frame_ca[idx]
        out.append((atom, xyz))
    return out


def write_fullatom_structures(path: str | Path, ca_records: Sequence[CARecord], start_ca: Array, frames: Sequence[Array], atom_records: Sequence[AtomRecord], keep_unmapped: bool = False, model_stride: int = 1) -> None:
    path=Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8') as out:
        midx=0
        serial_atoms=None
        for frame_idx, frame in enumerate(frames):
            if frame_idx % max(1, model_stride) != 0:
                continue
            midx += 1
            out.write(f"MODEL     {midx:4d}\n")
            decorated=decorate_frame_atoms(ca_records, start_ca, frame, atom_records, keep_unmapped=keep_unmapped)
            serial=1
            for atom, xyz in decorated:
                out.write(_format_atom_line(serial, atom, xyz)); serial += 1
            out.write('ENDMDL\n')
        out.write('END\n')


def write_selected_fullatom_frames(outdir: str | Path, ca_records: Sequence[CARecord], start_ca: Array, frames: Sequence[Array], atom_records: Sequence[AtomRecord], keep_unmapped: bool = False) -> None:
    from ba_arap_nma.output.writers import write_tsv
    outdir=Path(outdir)/'selected_fullatom_structures'; outdir.mkdir(parents=True, exist_ok=True)
    n=len(frames)
    idxs=sorted(set([0, max(0,n//4), max(0,n//2), max(0,3*n//4), max(0,n-1)]))
    rows=[]
    for idx in idxs:
        fn=f'frame_{idx:04d}_fullatom.pdb'
        write_fullatom_structures(outdir/fn, ca_records, start_ca, [frames[idx]], atom_records, keep_unmapped=keep_unmapped)
        rows.append({'frame':idx,'file':str(outdir/fn)})
    write_tsv(outdir/'selected_fullatom_structures.tsv', rows)


def has_fullatom_content(atom_records: Sequence[AtomRecord]) -> bool:
    return any(a.name.strip() != 'CA' for a in atom_records)
