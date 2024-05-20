#!/usr/bin/env python3

# Construct the supercell from unit cell

# supercell.py n1 n2 n3

# input: POSCAR
# output: POSCAR.new

from classes import POSCAR
import sys

if len(sys.argv)<4 :
    print("too few argv")
    print("expand by n1 n2 n3")
    sys.exit()
factor=[int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])]
factor_t=factor[0]*factor[1]*factor[2]

poscar1=POSCAR()
poscar1.fileread_vasp()

for i in range(3):
    for j in range(3):
        poscar1.lc[i][j]=poscar1.lc[i][j]*factor[i]
ia=0
apnew=[]
for i in range(poscar1.Ntype) :
    for j in range(poscar1.Naint[i]) :
        for ix in range(factor[0]) :
            for iy in range(factor[1]) :
                for iz in range(factor[2]) :
                    apnew.append([(poscar1.ap[ia][0]+ix)/factor[0],(poscar1.ap[ia][1]+iy)/factor[1],(poscar1.ap[ia][2]+iz)/factor[2]])
        ia=ia+1
poscar1.ap=apnew
for i in range(poscar1.Ntype) :
    poscar1.Naint[i]=poscar1.Naint[i]*factor_t
poscar1.Natom=poscar1.Natom*factor_t
poscar1.movetobox()

poscar1.filewrite_vasp()
