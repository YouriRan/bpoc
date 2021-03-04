import sys
import numpy as np
import MDAnalysis as mda

DUMPRESET = 20      # dump writing freq / xtc writing freq

u1 = mda.Universe(sys.argv[1], sys.argv[2])
u2 = mda.Universe(sys.argv[1], sys.argv[3])

# Determine last step to be read for first trajectory
maxts = len(u1.trajectory) - (len(u1.trajectory) % DUMPRESET)

out = sys.argv[1].rstrip('.data')

# Loop over both trajectories and write to .xtc
with mda.coordinates.XTC.XTCWriter(out + '_concat.xtc', n_atoms=u1.atoms.n_atoms) as f:
    for ts in u1.trajectory[:maxts]:
        f.write(u1.atoms)
    for ts in u2.trajectory[:]:
        f.write(u2.atoms)
