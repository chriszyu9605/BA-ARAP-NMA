from __future__ import annotations
import numpy as np
from ba_arap_nma.core import ca_core as core

def _norm_delta(start,target):
    d=(target-start).reshape(-1); n=np.linalg.norm(d)
    if n<1e-14: raise ValueError('The two states have zero coordinate difference')
    return d/n

def analyze_modal_space(start,target,segment_ids,max_modes=50,cutoff_A=15.0,summary_k=10):
    d=_norm_delta(start,target)
    linear=core.compute_low_modes(start,max_modes,cutoff_A)
    constrained=core.compute_constraint_compatible_modes(start,segment_ids,'C2',max_modes,cutoff_A)
    rows=[]; sl=0.0; sc=0.0
    for k in range(1,max_modes+1):
        sl+=float(linear.vectors[:,k-1]@d)**2; sc+=float(constrained.vectors[:,k-1]@d)**2
        rows.append({'k':k,'linear_cumulative_overlap':float(np.sqrt(sl)),'constrained_cumulative_overlap':float(np.sqrt(sc)),'delta_cumulative_overlap':float(np.sqrt(sc)-np.sqrt(sl))})
    rigid=core.rigid_body_basis(start)
    linear_proj=d-rigid@(rigid.T@d); eta_linear=float(linear_proj@linear_proj)
    j,_=core.build_local_distance_constraint_jacobian(start,segment_ids,'C2')
    combined=np.vstack([j,rigid.T]); P,rank=core._orthogonal_complement_basis(combined,3*len(start))
    eta=float(np.linalg.norm(P.T@d)**2)
    c_lin=rows[summary_k-1]['linear_cumulative_overlap']; c_con=rows[summary_k-1]['constrained_cumulative_overlap']
    raw=(target-start); global_rms=float(np.sqrt(np.mean(np.sum(raw*raw,axis=1))))
    scaled=raw/global_rms; s_loc=float(np.sqrt(np.mean((j@scaled.reshape(-1))**2)))
    summary={'n_ca':len(start),'segment_count':int(len(set(segment_ids.tolist()))),'local_constraint_count':int(j.shape[0]),'input_to_second_state_ca_rmsd_A':global_rms,'linear_space_eta':eta_linear,'constrained_space_eta':eta,'linear_space_maximum_overlap':float(np.sqrt(eta_linear)),'constrained_space_maximum_overlap':float(np.sqrt(eta)),'summary_k':summary_k,'linear_cumulative_overlap_at_summary_k':c_lin,'constrained_cumulative_overlap_at_summary_k':c_con,'delta_cumulative_overlap_at_summary_k':c_con-c_lin,'linear_leading_mode_concentration':c_lin*c_lin/eta_linear,'constrained_leading_mode_concentration':c_con*c_con/eta,'delta_leading_mode_concentration':c_con*c_con/eta-c_lin*c_lin/eta_linear,'local_distance_change_rms_per_1A_ca_rms_A':s_loc,'constraint_rank':rank-6,'constrained_nullity':P.shape[1]}
    return rows,summary,linear,constrained
