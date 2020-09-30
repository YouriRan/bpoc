import json
f = open("spica_top.json")
s= f.read()
d = json.loads(s)
f.close()

lst = []
for mol in d["topo"]:
    for t in d["topo"][mol]["type"]:
        if t not in [x[0] for x in lst]:
            c = d["topo"][mol]["type"].index(t)
            atomtype = []
            atomtype.append(d["topo"][mol]["type"][c])
            atomtype.append(d["topo"][mol]["charge"][c])
            atomtype.append(d["topo"][mol]["mass"][c])
            lst.append(atomtype)

#atp = open("~/bpoc/gromacs/SPICA.ff/atomtypes.atp")
atp = open("test")
for t in lst:
   print('{:>6s}{:>13f}'.format(t[0],t[2]))

atp.close()