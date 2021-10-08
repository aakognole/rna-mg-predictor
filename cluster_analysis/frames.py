# writes temporary pdb files for clustering purposes

import MDAnalysis
from MDAnalysis.analysis import align
from sys import argv

refpdb = ('../prep_system/system.nohyd.pdb')
dcd = str(argv[1])

u = MDAnalysis.Universe(refpdb, refpdb, dcd)
ref = MDAnalysis.Universe(refpdb)
align.alignto(u, ref, select="nucleic and name P", weights="mass")

full = u.select_atoms('all and not resname SOL')
for ts in u.trajectory[::100]:
    full.write('frame.'+str(ts.frame)+'.pdb')
