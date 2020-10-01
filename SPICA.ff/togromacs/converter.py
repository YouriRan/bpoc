# This script is written to transform the SPICA forcefield jsons to Gromacs
# parameter files.

# Written by:
# Youri A. Ran
# y.a.ran@umail.leidenuniv.nl

# The same structure as in MARTINI will be used.
# 2 ITP's and 1 TOP file are generated

### spica.itp ###
# Includes info on system, atomtypes and nonbond_params. 
#
# [defaults]
# ; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ
#
# [atomtypes]
# ; name mass charge ptype lj_1 lj_2
#
# [ nonbond_params ]
# ; i j funda lj_1 lj_2


### spica_mol.itp ###
# Includes molecule descriptors and constraint ifdefs
#
# [moleculetype]
# ; molname nrexcl
# nrexcl is the amount of bonds away at which moment to exclude neighbors
#
# [atoms]
# ; id type resnr residu atom cgnr charge
#
# [bonds]
# ; i j funct length force.c.
#
# [angles]
# ; i j k funct angle force.c.
#
# #ifdef ''
# [ position_restraints]
# ; i funct fcx fcy fcz


### system.top ###
# Includes .itp files and system info
#
# for every .itp file:
# #include "*.itp"
#
# [system]
# ; name
# 
# [molecules]
# ; name number

# Load SPICA Jsons
import json

f = open("spica_top.json")
s = f.read()
top = json.loads(s)
f.close()

f = open("spica_par.json")
s = f.read()
par = json.loads(s)
f.close()

# Unit conversion LAMMPS -> GROMACS
def unit(value, quantity):
    #                   LAMMPS          GROMACS         factor
    # length        r   A               nm              0.1
    # time          t   fs              ps              0.001
    # velocity      v   A/fs            nm ps-1         0.0001
    # mass          m   g/mol           u               1
    # temperature   T   K               K               1
    # pressure      P   atm.            bar             1.01325
    # energy        E,V kcal mol-1      kJ mol-1        4.2
    # force         F   kcal mol-1 A-1  kJ mol-1 nm-1   42
    # charge        q   e               e               1
    
    if quantity == 'r':
        return value*0.1
    elif quantity == 't':
        return value*0.001
    elif quantity == 'v':
        return value*0.0001
    elif quantity == 'm':
        return value
    elif quantity == 'T':
        return value
    elif quantity == 'P':
        return value*1.01325
    elif quantity == 'E' or quantity == 'V':
        return value*4.2
    elif quantity == 'F':
        return value*42
    elif quantity == 'q':
        return value
    else:
        print('Error: value not converted')
        exit(0)
    
def permutations(elements):
    if len(elements) <=1:
        yield elements
    else:
        for perm in permutations(elements[1:]):
            for i in range(len(elements)):
                # nb elements[0:1] works in both string and list contexts
                yield perm[:i] + elements[0:1] + perm[i:]

#############################################
#                 spica.itp                 #
#############################################

"""
TO DO
Convert epsilon and sigma to scaling constants for LJ (in the right units)
"""

defaults = [1, 1]
atomtypes = []
nonbond_params = []

# Extract all occurring cgtypes and their mass from spica_top
for mol in top['topo']:
    for cgtype in top['topo'][mol]['type']:
        if cgtype not in [x[0] for x in atomtypes]:
            index = top['topo'][mol]['type'].index(cgtype)
            l = []
            l.append(top['topo'][mol]['type'][index])
            l.append(top['topo'][mol]['mass'][index])
            atomtypes.append(l)
        
# Extract all interactions for pairs from spica_par
for interaction in par['params']:
    if interaction['param'] == 'pair':
        l = []
        l.append(interaction['types'])
        l.append(interaction['potential'])
        l.append(interaction['epsilon'])
        l.append(interaction['sigma'])
        nonbond_params.append(l)

# Write to file
spica_itp = open("./spica.itp",'w')

# Write defaults
spica_itp.write('[defaults]\n')
spica_itp.write('; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n')
spica_itp.write('{} {}\n\n'.format(defaults[0], defaults[1]))

# Write atomtypes
spica_itp.write('[atomtypes]\n')
spica_itp.write('; name mass charge ptype lj_1 lj_2\n')
for atomtype in atomtypes:
    # Charge and LJ are set to 0. 
    # Charge is specified in [moleculetype], LJ in [nonbond_params]
    spica_itp.write('{:>5s} {:>8.3f} 0 A 0 0\n'.format(atomtype[0],atomtype[1]))
spica_itp.write('\n')

# Write nonbond_params
spica_itp.write('[nonbond_params]\n')
spica_itp.write('; i j potential epsilon sigma\n')
spica_itp.write('; self terms\n')
for lj in nonbond_params:
    if lj[0][0] == lj[0][1]:
        spica_itp.write('{:>5s} {:>5s} {:>6s} {:>10.5f} {:>10.5f}\n'
            .format(lj[0][0], lj[0][1], lj[1], lj[2], lj[3]))
spica_itp.write('; cross terms\n')
for lj in nonbond_params:
    if lj[0][0] != lj[0][1]:
        spica_itp.write('{:>5s} {:>5s} {:>6s} {:>10.5f} {:>10.5f}\n'
            .format(lj[0][0], lj[0][1], lj[1], lj[2], lj[3]))
spica_itp.write('\n')

# Close files
spica_itp.close()

#############################################
#                 spica_mol.itp                 #
#############################################

# Uncomment this list and add molecules to only print those molecules
customlist = ['WAT', 'DOPC', 'POPC', 'DSPC']

moleculetypes = []
bonds = []
angles = []

spica_mol = open("spica_mol.itp", 'w')

# Extract molecule types from top
for molecule in top['topo']:
    # Add molecule to moleculetypes
    l = []
    l.append(molecule)
    l.append(top['topo'][molecule]['type'])
    l.append(top['topo'][molecule]['name'])
    l.append(top['topo'][molecule]['charge'])
    try:
        l.append(top['topo'][molecule]['bonds'])
    except:
        l.append('')
    try:
        l.append(top['topo'][molecule]['angles'])
    except:
        l.append('')
    moleculetypes.append(l)

# Extract all bonds and angles from par
for interaction in par['params']:
    if interaction['param'] == 'bond':
        l = []
        l.append(interaction['types'][0])
        l.append(interaction['types'][1])
        l.append(interaction['k'])
        l.append(interaction['r0'])
        bonds.append(l)
    elif interaction['param'] == 'angle':
        l = []
        l.append(interaction['types'][0])
        l.append(interaction['types'][1])
        l.append(interaction['types'][2])
        l.append(interaction['potential'])
        l.append(interaction['k'])
        l.append(interaction['theta0'])
        angles.append(l)

# Write to spica_mol.itp
for molecule in moleculetypes:
    if customlist:
        if molecule[0] not in customlist:
            continue    


    # Write [moleculetype]
    spica_mol.write('[moleculetype]\n')
    spica_mol.write('; molname nrexcl\n')
    spica_mol.write('{} 1\n\n'.format(molecule[0]))

    # Write [atoms]
    spica_mol.write('[atoms]\n')
    spica_mol.write('; id type resnr residu atom cgnr charge\n')
    cg_id = 1
    for cgtype in molecule[1]:
        spica_mol.write('{:>3d} {:>6s} {:>3d} {:>6s} {:>6s} {:>3d} {:>8.3f}\n'
            .format(cg_id, cgtype, 1, molecule[0], molecule[2][cg_id-1], cg_id, molecule[3][cg_id-1]))
        cg_id += 1
    spica_mol.write('\n')

    # Write [bonds]
    if molecule[4][0] != 'none':
        spica_mol.write('[bonds]\n')
        spica_mol.write('; i j funct length force.c.\n')
        
        # Loop over all bonds in molecule
        for bond in molecule[4]:
            i_t1 = molecule[2].index(bond[0]) + 1
            i_t2 = molecule[2].index(bond[1]) + 1
            t1 = molecule[1][i_t1-1]
            t2 = molecule[1][i_t2-1]

            check = False

            # Loop over all given harmonic bond potentials
            for b in bonds:
                for pair in permutations([t1,t2]):
                    if pair == [b[0],b[1]]:
                        bondtype = b
                        check = True

            # Uncomment to get list of unsupported bonds
            #if check == False:
                #print(f'Bond {t1} {t2} from molecule {molecule[0]} not found.')

            spica_mol.write('{:>3d} {:>3d} {:>3d} {:>8.3f} {:>8.3f}\n'
                .format(i_t1, i_t2, 1, unit(bondtype[3], 'r'), unit(bondtype[2],'F')))
    spica_mol.write('\n')

    # Write [angles]
    if molecule[5] == '' or molecule[5][0] == 'none':
        pass
    elif molecule[5][0] == 'auto':
        spica_mol.write('[angles]\n')
        spica_mol.write('i j k funct angle force.c.\n')

        # Loop over all bonds and check if types match
        for i in range(len(molecule[4])):
            for j in range(i+1, len(molecule[4])):
                check = False
                check2 = True
                bond1 = molecule[4][i]
                bond2 = molecule[4][j]
                angle = []
                #print(bond1, bond2)
                if bond1[0] == bond2[0]:
                    angle = [bond1[1], bond1[0], bond2[1]]
                elif bond1[1] == bond2[0]:
                    angle = [bond1[0], bond1[1], bond2[1]]
                elif bond1[0] == bond2[1]:
                    angle = [bond1[1], bond1[0], bond2[0]]
                elif bond1[1] == bond2[1]:
                    angle = [bond1[0], bond1[1], bond2[0]]
                else: 
                    check2 = False
                    
                if check2:
                    i_t1 = molecule[2].index(angle[0]) + 1
                    i_t2 = molecule[2].index(angle[1]) + 1
                    i_t3 = molecule[2].index(angle[2]) + 1
                    t1 = molecule[1][i_t1-1]
                    t2 = molecule[1][i_t2-1]
                    t3 = molecule[1][i_t3-1]

                    for a in angles:
                        # check angle permutations (fix middle atom)
                        if [t1, t2, t3] == [a[0], a[1], a[2]]:
                            angletype = a
                            check = True
                        elif [t3, t2, t1] == [a[0], a[1], a[2]]:
                            angletype = a
                            check = True

                    if check == True:
                        spica_mol.write('{:>3d} {:>3d} {:>3d} {:>8s} {:>8.3f} {:>8.3f}\n'
                            .format(i_t1, i_t2, i_t3, angletype[3], angletype[5], unit(angletype[4],'F')))
                    else:
                        # Uncomment to get list of unsupported angles
                        #print(f'Bond {t1} {t2} {t3} from molecule {molecule[0]} not found.')
                        spica_mol.write('{:>3d} {:>3d} {:>3d}\n'.format(i_t1, i_t2, i_t3))
    else:
        spica_mol.write('[angles]\n')
        spica_mol.write('; i j k funct angle force.c.\n')
        for angle in molecule[5]:
            i_t1 = molecule[2].index(angle[0]) + 1
            i_t2 = molecule[2].index(angle[1]) + 1
            i_t3 = molecule[2].index(angle[2]) + 1
            t1 = molecule[1][i_t1-1]
            t2 = molecule[1][i_t2-1]
            t3 = molecule[1][i_t3-1]

            check = False

            for a in angles:
                # check angle permutations (fix middle atom)
                if [t1, t2, t3] == [a[0], a[1], a[2]]:
                    angletype = a
                    check = True
                elif [t3, t2, t1] == [a[0], a[1], a[2]]:
                    angletype = a
                    check = True
            
            if check == True:
                spica_mol.write('{:>3d} {:>3d} {:>3d} {:>8s} {:>.3f} {:>.3f}\n'
                    .format(i_t1, i_t2, i_t3, angletype[3], angletype[5], unit(angletype[4],'F')))
            else:
                # Uncomment to get list of unsupported angles
                #print(f'Bond {t1} {t2} {t3} from molecule {molecule[0]} not found.')
                pass

    spica_mol.write('\n')
    spica_mol.write(';----------------------------------------------------------------\n\n')
spica_mol.close()

f = open("angle.txt",'w')
for angle in angles:
    f.write('{} {} {} {} {} {}\n'.format(angle[0], angle[1], angle[2], angle[3], angle[4], angle[5]))
f.close()


### spica_lip.itp ###
# Includes molecule descriptors and constraint ifdefs
#
# [moleculetype]
# ; molname nrexcl
# nrexcl is the amount of bonds away at which moment to exclude neighbors
#
# [atoms]
# ; id type resnr residu atom cgnr charge
#
# [bonds]
# ; i j funct length force.c.
#
# [angles]
# ; i j k funct angle force.c.
#
# #ifdef ''
# [ position_restraints]
# ; i funct fcx fcy fcz