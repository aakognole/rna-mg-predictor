import MDAnalysis as mda

u = mda.Universe('system.pdb')
nohyd = u.select_atoms('all and not name H*')
nohyd.write('system.nohyd.pdb')

rna = u.select_atoms('nucleic and not name H*')
n = rna.n_atoms
f1 = open('natoms.txt','w')
f1.write(str(n))
f1.close()

rna = u.select_atoms('nucleic')
com = rna.center_of_mass()
f2 = open('rna_com.txt','w')
f2.write("%.1f,%.1f,%.1f"%(com[0]/10,com[1]/10,com[2]/10))
f2.close()
