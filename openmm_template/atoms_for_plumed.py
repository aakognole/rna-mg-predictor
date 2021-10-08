# Plumed scripts are based on atom numbers
# Here, MDAnalysis is used to get the set of atoms to define CV

from __future__ import print_function
import MDAnalysis as mda

print('Note: Maybe delete CRYST1 line from pdb')
print('Note: Do not copy the last comma to plumed script')
print()

pdb = ('../prep_system/system.pdb')
u = mda.Universe(pdb)

g1 = u.select_atoms("nucleic and (resid 1:22) and not name H*")
print('ATOMS', end='=')
for i in g1.atoms.indices:
    print(i, end=',') 

print()
print()

g2 = u.select_atoms("nucleic and (resid 23:34) and not name H*")
print('ATOMS', end='=')
for i in g2.atoms.indices:
    print(i, end=',')
