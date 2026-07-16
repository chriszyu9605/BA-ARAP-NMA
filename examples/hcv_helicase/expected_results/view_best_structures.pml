reinitialize
load paired_second_state_ca.pdb, target
load ba_arap_nma_best_ca.pdb, ba
load linear_ca_anm_best_ca.pdb, linear
hide everything, all
set cartoon_trace_atoms, on
show cartoon, target or ba or linear
show sticks, target or ba or linear
color green, target
color orange, ba
color marine, linear
orient target or ba or linear
bg_color white
