# To reduce the storage, this file removes hydrogens from the trajectory
# Also, the last coordinate and trajectory is aligned to the initial
# structure and rmsd is calculated.

import MDAnalysis
import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

pdb = ('solution.<run>.prod.<curr>.pdb')
dcd = ('solution.<run>.prod.<curr>.dcd')
nohpdb = ('../prep_system/system.nohyd.pdb')
ref = MDAnalysis.Universe(nohpdb)

u = MDAnalysis.Universe(pdb)
nohyd = u.select_atoms("all and not name H*")
nohyd.write('temp.nohyd.pdb')

u = MDAnalysis.Universe('temp.nohyd.pdb')
align.alignto(u,ref,select="nucleic and not name H*",weights="mass")
nohyd = u.select_atoms("all and not name H*")
nohyd.write('solution.<run>.prod.<curr>.nohyd.pdb')

u = MDAnalysis.Universe(pdb, dcd)
nohyd = u.select_atoms("all and not name H*")
with MDAnalysis.Writer("temp.nohyd.dcd", nohyd.n_atoms) as W:
    for ts in u.trajectory:
        W.write(nohyd)

noh = MDAnalysis.Universe(nohpdb, 'temp.nohyd.dcd')
alignment = align.AlignTraj(noh, ref, select="nucleic and not name H*", filename='solution.<run>.prod.<curr>.nohyd.dcd')
alignment.run()

a = MDAnalysis.Universe(nohpdb, nohpdb, 'solution.<run>.prod.<curr>.nohyd.dcd')
rna = a.select_atoms("nucleic and not name H*")
refrna = ref.select_atoms("nucleic and not name H*")
RMSD = []
for ts in a.trajectory:
    R = rmsd(rna.positions, refrna.positions)
    RMSD.append((ts.frame, R))
RMSD = np.array(RMSD)
np.savetxt('rmsd.<run>.<curr>.txt', RMSD, fmt='%s')
