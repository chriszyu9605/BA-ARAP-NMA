from __future__ import annotations
from pathlib import Path
from typing import Sequence
import csv,json
import numpy as np
from numpy.typing import NDArray
from ba_arap_nma.core.ca_core import CARecord
from .fullatom import AtomRecord,write_fullatom_structures,write_selected_fullatom_frames,has_fullatom_content
Array=NDArray[np.float64]

def write_tsv(path,rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if not rows: path.write_text('',encoding='utf-8'); return
    fields=[]
    for row in rows:
        for k in row:
            if k not in fields: fields.append(k)
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

def write_json(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=True),encoding='utf-8')

def _atom_line(serial,rec,xyz):
    chain=rec.chain if rec.chain!='_' else ' '
    return f"ATOM  {serial:5d}  CA  {rec.resname:>3s} {chain:1s}{rec.resseq:4d}{rec.icode:1s}   {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00           C\n"

def _bonds(records):
    return [(i+1,i+2) for i in range(len(records)-1) if records[i].chain==records[i+1].chain and records[i+1].resseq-records[i].resseq==1]

def write_ca_structures(path,records,frames,conect=False):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True); bonds=_bonds(records)
    with path.open('w',encoding='utf-8') as out:
        for model,coords in enumerate(frames,1):
            out.write(f'MODEL     {model:4d}\n')
            for serial,(rec,xyz) in enumerate(zip(records,coords),1): out.write(_atom_line(serial,rec,xyz))
            if conect:
                for i,j in bonds: out.write(f'CONECT{i:5d}{j:5d}\n')
            out.write('ENDMDL\n')
        out.write('END\n')

def write_single(path,records,coords):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8') as out:
        for serial,(rec,xyz) in enumerate(zip(records,coords),1): out.write(_atom_line(serial,rec,xyz))
        for i,j in _bonds(records): out.write(f'CONECT{i:5d}{j:5d}\n')
        out.write('END\n')

def write_ca_pymol(outdir):
    write(Path(outdir)/'view_ca_structures.pml',"reinitialize\nload ca_trace_structures.pdb, structures\nhide everything, all\nset all_states, on\nset cartoon_trace_atoms, on\nshow cartoon, structures\nshow sticks, structures\nset stick_radius, 0.18\nset cartoon_tube_radius, 0.45\norient structures\nbg_color white\n# Use mplay to animate states.\n")

def write(path,text): Path(path).write_text(text,encoding='utf-8')

def selected_indices(n): return sorted(set([0,n//4,n//2,3*n//4,n-1]))

def write_branch_outputs(outdir,records,frames,rows,summary,config,fullatom_records=None):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    write_ca_structures(outdir/'ca_structures.pdb',records,frames,False)
    write_ca_structures(outdir/'ca_trace_structures.pdb',records,frames,True)
    write_ca_pymol(outdir)
    sd=outdir/'selected_structures'; sd.mkdir(exist_ok=True)
    selected_rows=[]
    for idx in selected_indices(len(frames)):
        fn=f'frame_{idx:04d}.pdb'; write_single(sd/fn,records,frames[idx]); selected_rows.append({'frame':idx,'file':fn})
    write_tsv(sd/'selected_structures.tsv',selected_rows)
    full_status='not_generated'
    if fullatom_records and has_fullatom_content(fullatom_records):
        start=np.stack([r.coord for r in records])
        write_fullatom_structures(outdir/'fullatom_decorated_structures.pdb',records,start,frames,fullatom_records)
        write_selected_fullatom_frames(outdir,records,start,frames,fullatom_records)
        write(outdir/'view_fullatom_decorated.pml',"reinitialize\nload fullatom_decorated_structures.pdb, structures\nhide everything, all\nshow cartoon, structures\nshow sticks, (structures and name CA+N+C+O)\nset all_states, on\norient structures\nbg_color white\n# Use mplay to animate states.\n")
        full_status='generated_from_input_pdb_atom_records'
    summary=dict(summary); summary['fullatom_decorated_output']=full_status
    config=dict(config); config['fullatom_decorated_output']=full_status
    write_tsv(outdir/'frame_metrics.tsv',rows); write_json(outdir/'branch_summary.json',summary); write_json(outdir/'run_config.json',config)
