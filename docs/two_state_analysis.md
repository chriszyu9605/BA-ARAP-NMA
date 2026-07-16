# Two-state analysis

The `two-state` command maps common protein C-alpha residues between the two PDB files, aligns the paired second state to the input, calculates the standard and constrained modal bases, and writes cumulative-overlap and constraint-space quantities. The paired second state is not used to generate BA-ARAP-NMA coordinates.

Theory outputs include cumulative overlap for K=1 to the requested analysis-mode count, complete constrained-space displacement retention (eta), maximum attainable overlap (sqrt(eta)), leading-mode concentration (rho_K), and the first-order local-distance change per 1 A global C-alpha RMS displacement (s_loc).
