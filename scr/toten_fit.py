#!/bin/python3

# print the distortion coordinate (out-of-plane displacement of the central atom)
# and the total energy
from vasp import filein

N=41

# True: total energy
# False: band energy
is_toten=False

poscar_0=filein.POSCAR("POSCAR_i")


energy=[]
pos=[]

for i in range(N) :

    if is_toten :
        f1=open(str(i+1)+"/OSZICAR")
        line=f1.readlines()
        f1.close()
        word=line[-1].split()
        energy.append(float(word[2]))
    else :
        eigenval_1=filein.EIGENVAL(str(i+1)+"/EIGENVAL")
        energy0=[]
        for ispin in range(2) :
            for ib in range(140,145) :
                energy0.append(eigenval_1.eig[0][ib][ispin])
        energy.append(energy0)
        del energy0

    poscar_1=filein.POSCAR(str(i+1)+"/POSCAR")
    

##  Q = Renormalization factor:
#    pos.append(poscar_1.total_distance(poscar_0))

#   Q = displacement of the central factor
    pos.append(poscar_1.lc[2][2]*poscar_1.ap[35][2])

    del poscar_1

if is_toten :
    energy_zero=energy[0]
    for i in range(N):
        energy[i]=energy[i]-energy_zero
else :
    energy_zero=energy[0][2]
    for i in range(N):
        for j in range(len(energy[i])) :
            energy[i][j]=energy[i][j]-energy_zero

energy_out=[]
pos_out=[]
for i in range(N-1):
    pos_out.append(-pos[-i-1])
    energy_out.append(energy[-i-1])

for i in range(N):
    pos_out.append(pos[i])
    energy_out.append(energy[i])

if is_toten :
    for i in range(2*N-1):
        print(pos_out[i],energy_out[i])
else :
    fout=open("post.out","w")
    for i in range(2*N-1):
        fout.write("%.6f"%pos_out[i])
        for j in range(len(energy_out[i])) :
            fout.write(" %.6f"%energy_out[i][j])
        fout.write("\n")
