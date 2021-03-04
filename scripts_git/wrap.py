import sys
import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis import align
import MDAnalysis.transformations as trans
from MDAnalysis.analysis.leaflet import LeafletFinder

u = mda.Universe(sys.argv[1] + '.data', sys.argv[1] + '.xtc')
ref = mda.Universe(sys.argv[1] + '.data')

# Get box dimensions from first frame
b = u.trajectory[0].dimensions[0:3]
bbox = np.array([0, 0, 0, b[0], b[1], b[2]]).reshape(2,3)

# Align all frames to first frame on PH beads
alignment = mda.analysis.align.AlignTraj(u, ref, select='type 2')
alignment.run()

# Translate box to centre buckle
translator = [0.5*bbox[1, 0], 0.5*bbox[1, 1], 0.5*bbox[1, 2]]

# Feed transformations to trajectory
wf = [trans.translate(translator), trans.wrap(u.atoms, compound='residues')]
u.trajectory.add_transformations(*wf)

# Write output
u.atoms.write(sys.argv[1] + '_wrap.xtc', frames=u.trajectory[:])
