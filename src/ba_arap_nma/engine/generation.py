from __future__ import annotations
from dataclasses import dataclass,asdict
import time
import numpy as np
from numpy.typing import NDArray
from ba_arap_nma.core import ca_core as core
from ba_arap_nma.core.arap import arap_tangent_preserving_project
from ba_arap_nma.core.structure import Structure
from ba_arap_nma.analysis.geometry import local_geometry_metrics
Array=NDArray[np.float64]

@dataclass
class GenerationConfig:
    mode_rank:int=1; orientation:int=1; cutoff_A:float=15.0; max_amplitude_A:float=20.0; amplitude_step_A:float=0.25

@dataclass
class BranchResult:
    branch_id:str; frames:list[Array]; rows:list[dict[str,object]]; summary:dict[str,object]; config:dict[str,object]

def normalize_mode(v,n):
    x=np.asarray(v,float).reshape(n,3); rms=float(np.sqrt(np.mean(np.sum(x*x,axis=1))))
    if rms<1e-14: raise ValueError('Mode has zero C-alpha RMS amplitude')
    return x/rms

def _target_metrics(initial,coords,target):
    if target is None: return {'target_rmsd_A':float('nan'),'transition_coverage':float('nan')}
    r0=core.ca_rmsd(initial,target,True); r=core.ca_rmsd(coords,target,True)
    return {'target_rmsd_A':r,'transition_coverage':(r0-r)/r0 if r0>1e-12 else 1.0}

def generate_branch(structure:Structure,cfg:GenerationConfig,mode_vector:Array,method:str='ba_arap_nma',target:Array|None=None):
    t0=time.time(); initial=structure.coords.copy(); step=cfg.amplitude_step_A
    if step<=0 or cfg.max_amplitude_A<0: raise ValueError('Amplitude values must be non-negative and step must be positive')
    count=int(round(cfg.max_amplitude_A/step))
    if abs(count*step-cfg.max_amplitude_A)>1e-8: raise ValueError('max amplitude must be an integer multiple of amplitude step')
    mode=normalize_mode(np.asarray(mode_vector)*float(cfg.orientation),len(initial))
    branch_id=f"{method}_mode{cfg.mode_rank:03d}_{'plus' if cfg.orientation>0 else 'minus'}"
    frames=[]; rows=[]; prev=initial.copy()
    for frame in range(count+1):
        amp=frame*step
        if frame==0: coords=initial.copy()
        else:
            predictor=initial+amp*mode
            if method=='ba_arap_nma':
                coords,_=arap_tangent_preserving_project(initial,predictor,structure.segment_ids,component_gauge_target='predictor')
                coords,_,_=core.kabsch_align(coords,prev)
            elif method=='linear_ca_anm': coords=predictor
            else: raise ValueError(method)
        geom=local_geometry_metrics(initial,coords,structure.segment_ids)
        row={'branch_id':branch_id,'method':method,'mode_rank':cfg.mode_rank,'orientation':cfg.orientation,'frame':frame,'requested_amplitude_A':amp,'ca_rmsd_from_input_A':core.ca_rmsd(coords,initial,True)}
        row.update(geom); row.update(_target_metrics(initial,coords,target)); rows.append(row); frames.append(coords.copy()); prev=coords
    finite=[r for r in rows if np.isfinite(float(r['target_rmsd_A']))]
    best=min(finite,key=lambda r:float(r['target_rmsd_A'])) if finite else None
    summary={'branch_id':branch_id,'method':method,'mode_rank':cfg.mode_rank,'orientation':cfg.orientation,'frame_count':len(frames),'runtime_s':time.time()-t0,'max_requested_amplitude_A':cfg.max_amplitude_A,'amplitude_step_A':step}
    if best:
        summary.update({'best_target_frame':best['frame'],'best_target_requested_amplitude_A':best['requested_amplitude_A'],'best_target_rmsd_A':best['target_rmsd_A'],'best_transition_coverage':best['transition_coverage'],'best_input_relative_virtual_angle_mae_deg':best['virtual_angle_mae_deg']})
    return BranchResult(branch_id,frames,rows,summary,asdict(cfg)|{'method':method,'fixed_initial_mode':True,'recompute_modes_each_frame':False})
