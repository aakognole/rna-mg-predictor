import MDAnalysis as mda
u = mda.Universe('system.pdb')
nohyd = u.select_atoms('all and not name H*')
nohyd.write('system.nohyd.pdb')
rna = u.select_atoms('nucleic and not name H*')
n = rna.n_atoms
f = open('natoms.txt','w')
f.write(str(n))
f.close()
