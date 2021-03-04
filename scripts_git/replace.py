import numpy as np
import random
import sys

pdb_in = open("4j_3.pdb", 'r')
phi_dogl = 0.15
pdb_out = open(str(int(phi_dogl * 100)) + ".pdb", 'w')
sections = [500,100,100,100,100,100,100,100,100]

luckynumbers = []
count = 1

for sect in sections:
    luckynumbers.append(random.sample(range(count, count+sect),int(phi_dogl*sect)))
    count += sect    

# flatten list
l = []
for i in luckynumbers:
    for j in i:
        l.append(j)
luckynumbers = l

def pdbreader(line):
    columns = []
    columns.append(str(line[0:6]).strip())
    columns.append(str(line[6:11]).strip())
    columns.append(str(line[12:16]).strip())
    columns.append(str(line[16:21]).strip())
    columns.append(str(line[21:22]).strip())
    columns.append(str(line[22:26]).strip())
    columns.append(str(line[30:38]).strip())
    columns.append(str(line[38:46]).strip())
    columns.append(str(line[46:54]).strip())
    columns.append(str(line[54:60]).strip())
    columns.append(str(line[60:66]).strip())
    columns.append(str(line[67:77]).strip())
    return columns

# New .pdb is written in 3 steps
# First read and write all lines containing headers and DOPC not in luckynumbers
# Second read and write all lines containing DOPC in luckynumbers, skip NC3, replace PH with OG
# Third read and write all lines containing water and END
# n keeps track of new atom ID number
n = 1
for line in pdb_in:
    columns = pdbreader(line)
    if columns[0] == 'ATOM' and columns[3] == 'DOPC' and int(columns[5]) not in luckynumbers:
        pdb_out.write('{:<5s}{:>6d} {:<4s}{:<5s}{:<1s}{:>4s}    {:>8s}{:>8s}{:>8s}{:>6s}{:>6s}    {:>6s}\n'
            .format(columns[0], n, columns[2], columns[3], columns[4], columns[5], columns[6], columns[7], columns[8], columns[9], columns[10],columns[11]))
        n += 1
    elif columns[0] in ['TITLE','HEADER','REMARK']:
        pdb_out.write(line)

pdb_in.seek(0)

for line in pdb_in:
    columns = pdbreader(line)
    if columns[0] == 'ATOM' and columns[3] == 'DOPC' and int(columns[5]) in luckynumbers:
        if columns[2] in ['GL', 'EST1', 'EST2', 'C11','C12','C13','C14','C15','C16','C21','C22','C23','C24','C25','C26']:
            pdb_out.write('{:<5s}{:>6d} {:<4s}{:<5s}{:<1s}{:>4s}    {:>8s}{:>8s}{:>8s}{:>6s}{:>6s}    {:>6s}\n'
                .format(columns[0], n, columns[2], 'DOGL', columns[4], columns[5], columns[6], columns[7], columns[8], columns[9], columns[10],columns[11]))
            n += 1
        elif columns[2] == 'PH':
            pdb_out.write('{:<5s}{:>6d} {:<4s}{:<5s}{:<1s}{:>4s}    {:>8s}{:>8s}{:>8s}{:>6s}{:>6s}    {:>6s}\n'
                .format(columns[0], n, 'OG', 'DOGL', columns[4], columns[5], columns[6], columns[7], columns[8], columns[9], columns[10],columns[11]))
            n += 1

pdb_in.seek(0)

for line in pdb_in:
    columns = pdbreader(line)
    if columns[0] == 'ATOM' and columns[3] == 'TIP3':
        pdb_out.write('{:<5s}{:>6d} {:<4s}{:<5s}{:<1s}{:>4s}    {:>8s}{:>8s}{:>8s}{:>6s}{:>6s}    {:>6s}\n'
            .format(columns[0], n, columns[2], columns[3], columns[4], columns[5], columns[6], columns[7], columns[8], columns[9], columns[10],columns[11]))
        n += 1
    elif columns[0] == 'END':
        pdb_out.write(line)
