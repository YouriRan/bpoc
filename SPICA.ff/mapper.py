#!/usr/bin/python3
# SPICA mapper aa.pdb -> cg.pdb
# Written by:
# Youri Ran
# y.a.ran@umail.leidenuniv.nl

# Input is read from stdin
import fileinput
pdb_in = fileinput.input()

# Output is printed (stdout) so redirect to .pdb file

import json
import os
dirname = os.path.dirname(__file__)
spica_top = os.path.join(dirname, 'spica/spica_top_dog.json')
f = open(spica_top, 'r')
s = f.read()
top = json.loads(s)
dopc = top['topo']['DOPC']
wat = top['topo']['WAT']
dogl = top['topo']['DOGL']
cla = top['topo']['CLA']
sod = top['topo']['SOD']
f.close()

# Save all molecules as a list of atoms in list molecules
molecules = []
atoms = []
molcounter = 0

for line in pdb_in:
    columns = line.split()
    if columns[0] == 'REMARK':
        for word in columns:
            print(word, end=' ')
        print()
    elif columns[0] != 'ATOM':
        continue
    else:
        # Check if new molecule
        if int(columns[4]) != molcounter:
            molecules.append(atoms)
            molcounter = int(columns[4])
            atoms = []
        
        atoms.append(columns)
molecules.append(atoms)

# print own remark
print("REMARK MAPPED TO SPICA COARSE-GRAINED FORCE FIELD BY YOURI RAN")

# For all molecules print all atoms and coordinates
molcounter = 1
beadcounter = 1
for mol in molecules[1:]:
    molname = mol[0][3]

    # Check which molecule and load map and names of cgbeads
    if molname == 'DOPC':
        beadlist = [[] for x in range(len(dopc['map']))]
        cgmaps = dopc['map']
        names = dopc['name']
        moltype = 'MEMB'
    elif molname == 'TIP3':
        molname = 'WAT'
        beadlist = [[] for x in range(len(wat['map']))]
        cgmaps = wat['map']
        names = wat['name']
        moltype = 'WAT'
    elif molname == 'DOGL':
        beadlist = [[] for x in range(len(dogl['map']))]
        cgmaps = dogl['map']
        names = dogl['name']
        moltype = 'MEMB'
        

    # make list with atoms and positions per cg particle
    for atom in mol:
        for cgmap in cgmaps:
            if atom[2] in cgmap:
                index = cgmaps.index(cgmap)
                beadlist[index].append([atom[2], atom[5], atom[6], atom[7]])

    cgbeads = []
    # Map to cg particles by centre-of-mass
    for bead in beadlist:
        x, y, z = 0, 0, 0
        mass = 0
        for atom in bead:
            # Determine mass of atom (all tip3 waters weigh 18)
            if molname == 'DOPC' or molname == 'DOGL':
                if atom[0][0] == 'C':
                    m = 12
                elif atom[0][0] == 'H':
                    m = 1
                elif atom[0][0] == 'N':
                    m = 14
                elif atom[0][0] == 'P':
                    m = 31
                elif atom[0][0] == 'O':
                    m = 16
            elif molname == 'WAT':
                m = 18

            # Add mass-weighted coordinate
            x += m*float(atom[1])
            y += m*float(atom[2])
            z += m*float(atom[3])
            mass += m

        # Normalize by total mass
        x = x/mass
        y = y/mass
        z = z/mass

        # Add cgbead index and coordinate
        index = beadlist.index(bead)
        cgbeads.append([index, x, y, z])

    # Write all beads
    for bead in cgbeads:
        print('{:<5s}{:>6d} {:<4s}{:>5s}{:>5d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}    {:>6s}'
            .format('ATOM', beadcounter, names[bead[0]], molname, molcounter, bead[1], bead[2], bead[3], 1.00, 0.00, moltype))
        beadcounter += 1
    molcounter +=1

print('{:<5s}{:>6d} {:<4s}{:>5s}{:>5d}'
    .format('TER', beadcounter, '', molname, molcounter))
print('{:>5s}'.format('END'))


