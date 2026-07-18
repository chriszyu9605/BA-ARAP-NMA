from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
Array=NDArray[np.float64]

def _pairs(segment_ids, gap):
    return [(i,i+gap) for i in range(len(segment_ids)-gap) if segment_ids[i]==segment_ids[i+gap]]

def _angle(a,b,c):
    u=a-b; v=c-b; nu=np.linalg.norm(u); nv=np.linalg.norm(v)
    if nu<1e-12 or nv<1e-12: return np.nan
    return np.degrees(np.arccos(np.clip(np.dot(u,v)/(nu*nv),-1.0,1.0)))

def local_geometry_metrics(reference: Array, coords: Array, segment_ids):
    edges=_pairs(segment_ids,1)
    d0=np.array([np.linalg.norm(reference[i]-reference[j]) for i,j in edges])
    d=np.array([np.linalg.norm(coords[i]-coords[j]) for i,j in edges])
    dev=d-d0
    triples=[(i,i+1,i+2) for i in range(len(segment_ids)-2) if segment_ids[i]==segment_ids[i+1]==segment_ids[i+2]]
    a0=np.array([_angle(reference[i],reference[j],reference[k]) for i,j,k in triples])
    a=np.array([_angle(coords[i],coords[j],coords[k]) for i,j,k in triples])
    adev=np.abs(a-a0); adev=adev[np.isfinite(adev)]
    return {
      'mean_squared_adjacent_ca_distance_deviation_A2': float(np.mean(dev*dev)) if len(dev) else float('nan'),
      'adjacent_ca_distance_mae_A': float(np.mean(np.abs(dev))) if len(dev) else float('nan'),
      'adjacent_distance_outliers_gt_0p5_A': int(np.sum(np.abs(dev)>0.5)),
      'adjacent_distance_outliers_gt_1p0_A': int(np.sum(np.abs(dev)>1.0)),
      'virtual_angle_mae_deg': float(np.mean(adev)) if len(adev) else float('nan'),
      'virtual_angle_outliers_gt_10_deg': int(np.sum(adev>10.0)) if len(adev) else 0,
    }
