import json

f = open("spica_top.json")
s = f.read()
top = json.loads(s)
f.close()

dspc = top['topo']['DSPC']

c36 = open("/home/youriran/Downloads/dspc_min.xyz", 'r')

beadlist = [[] for x in range(len(dspc['map']))]

for line in c36:
    l = line.split()
    if len(l) == 4:
        check = False
        for m in dspc['map']:
            if l[0] in m:
                index = dspc['map'].index(m)
                beadlist[index].append([l[1],l[2],l[3]])

spica = open("spica_dspc.xyz",'w')
spica.write(str(len(beadlist)))
spica.write('\nDSPC\n')    

for bead in beadlist:
    x, y, z = 0, 0, 0
    for coord in bead:
        x += float(coord[0])
        y += float(coord[1])
        z += float(coord[2])
    x = x / len(bead)
    y = y / len(bead)
    z = z / len(bead)
    
    t = dspc['name'][beadlist.index(bead)]

    spica.write('{:>6s} {:>6.3f} {:>6.3f} {:>6.3f}\n'.format(t, x, y, z))

spica.close()

a = 1.0
f = open("w.xyz",'w')
f.write('1\nW\n')
f.write('{:>6s} {:>6.3f} {:>6.3f} {:>6.3f}\n'.format('W', a,a,a))
f.close()