import math

f1=open("POSCAR","r")
line=f1.readlines()
f1.close()

lc=[[],[],[]]
lf=float(line[1].split()[0])
for i in range(3) :
    word0=line[i+2].split()
    for j in range(3) :
        lc[i].append(float(word0[j])*lf)
satom=line[5].split()
word0=line[6].split()
Natom=0
iatom=[]
atom=[]
for i in range(len(word0)) :
    iatom.append(int(word0[i]))
    Natom=Natom+iatom[i]
for i in range(Natom) :
    atom1=[]
    word0=line[8+i].split()
    for j in range(3) :
        atom1.append(float(word0[j]))
        if atom1[j]==-0.0 :
            atom1[j]=0.0
    atom.append(atom1)

kgrid=[]
for i in range(3) :
    kgrid.append(math.floor(30.0/(lc[i][0]**2+lc[i][1]**2+lc[i][2]**2)**0.5))

f2=open("qe.st","w")
f2.write("CELL_PARAMETERS angstrom\n")
for i in range(3) :
    f2.write(f"  {lc[i][0]:.12f}  {lc[i][1]:.12f}  {lc[i][2]:.12f}\n")
f2.write("ATOMIC_SPECIES\n")
for i in range(len(satom)) :
    f2.write("  "+satom[i]+"  1.0  PP\n")
f2.write("ATOMIC_POSITIONS crystal\n")
ij=0
ik=0
for i in range(len(atom)) :
    f2.write(f"  {satom[ij]:2s}  {atom[i][0]:.16f}  {atom[i][1]:.16f}  {atom[i][2]:.16f}\n")
    ik=ik+1
    if ik==iatom[ij] :
        ij=ij+1
        ik=0
f2.write("K_POINTS automatic\n")
f2.write(f"  {kgrid[0]:d}  {kgrid[1]:d}  {kgrid[2]:d} 0 0 0\n")
f2.close()
