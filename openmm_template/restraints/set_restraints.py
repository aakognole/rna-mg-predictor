# This script writes a file to be used add restraints by CHARMM-GUI OpenMM inputs
import MDAnalysis as mda

u = mda.Universe('../../prep_system/system.pdb')

sel = u.select_atoms("(protein or nucleic) and not name H* *' O*P")
BB = u.select_atoms('backbone or nucleicbackbone')
SC = []
for n in sel:
    if n not in BB:
        SC.append(n)

f = open('prot_pos.txt','w')
rna = u.select_atoms("(protein or nucleic) and not name H*")
for i in rna.atoms:
    if i in BB:
        #print("%d BB" % (i.id))
        f.write("%d BB \n" % i.id)
for i in rna.atoms:
    if i in SC:
        #print("%d SC" % (i.id))
        f.write("%d SC \n" % i.id)
f.close()
