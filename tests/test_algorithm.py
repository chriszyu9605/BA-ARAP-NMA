import numpy as np
from pathlib import Path
from ba_arap_nma.core.structure import load_structure
from ba_arap_nma.core import ca_core as core
from ba_arap_nma.engine.generation import GenerationConfig,generate_branch
ROOT=Path(__file__).parent

def test_constraints_and_rigid_body():
    st=load_structure(ROOT/'data/toy_input.pdb')
    m=core.compute_constraint_compatible_modes(st.coords,st.segment_ids,'C2',1,15.0)
    assert m.rigid_projection_max_abs<1e-8
    assert m.constraint_projection_max_abs<1e-8

def test_81_frames_and_target_independence():
    st=load_structure(ROOT/'data/toy_input.pdb')
    basis=core.compute_constraint_compatible_modes(st.coords,st.segment_ids,'C2',1,15.0)
    cfg=GenerationConfig(1,1,15.0,20.0,0.25)
    a=generate_branch(st,cfg,basis.vectors[:,0],'ba_arap_nma',None)
    target=st.coords+0.1
    b=generate_branch(st,cfg,basis.vectors[:,0],'ba_arap_nma',target)
    assert len(a.frames)==81
    for x,y in zip(a.frames,b.frames): assert np.allclose(x,y,atol=1e-12)
