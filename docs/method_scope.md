# Method scope

The implementation uses a 15 Å single-parameter C-alpha anisotropic network by default, analytically removes three translations and three rotations, constrains first-order changes in adjacent and next-nearest C-alpha distances, re-diagonalizes the Hessian in the joint admissible space, and applies a start-referenced co-rotational ARAP correction at each requested amplitude. Modes are not recomputed frame by frame.
