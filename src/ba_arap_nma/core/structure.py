from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from numpy.typing import NDArray
from . import ca_core as core

Array = NDArray[np.float64]

@dataclass
class Structure:
    records: list[core.CARecord]
    coords: Array
    segment_ids: NDArray[np.int64]
    parse_audit: core.ParseAudit

def segment_ids(records: list[core.CARecord]) -> NDArray[np.int64]:
    seg=np.zeros(len(records),dtype=np.int64); cur=0
    for i in range(1,len(records)):
        a,b=records[i-1],records[i]
        if b.chain!=a.chain or (b.resseq-a.resseq)>1 or (b.resseq-a.resseq)<0:
            cur+=1
        seg[i]=cur
    return seg

def load_structure(pdb_path: str|Path, chain: str|None=None) -> Structure:
    records,audit=core.parse_ca_records_with_audit(pdb_path,chain)
    coords=np.stack([r.coord for r in records])
    if len(records)<4:
        raise ValueError('At least four protein C-alpha atoms are required')
    return Structure(records,coords,segment_ids(records),audit)
