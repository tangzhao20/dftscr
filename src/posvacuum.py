#!/usr/bin/env python3
#
# This script adds the vacuum layer to the structure.
#
# python posvacuum.py z_vac(A)
# where z_vac is in Angstrom.
#
# Input: POSCAR

from classes import POSCAR
import sys

if len(sys.argv)<2 :
    print("Error: too few arguments")
    print("python posvacuum.py z_vac(A)")
    sys.exit()
else :
    z_vac=float(sys.argv[1])

if z_vac<0 :
    print("Error: z_vac < 0")
    sys.exit()

poscar1=POSCAR()
poscar1.fileread_vasp()

if max(abs(poscar1.lc[2][0]), abs(poscar1.lc[2][1]))>0.000001 :
    print("Error: a3 doesn't parllel to z")
    sys.exit()
a3_new=poscar1.lc[2][2]+z_vac
factor=(a3_new/poscar1.lc[2][2])

poscar1.movetobox()

ia=0
newatoms=[]
for itype in range(poscar1.Ntype) :
    newatoms0=[]
    for iaint in range(poscar1.Naint[itype]) :
        if poscar1.ap[ia][2]<0.00000001 :
            newatoms0.append([poscar1.ap[ia][0],poscar1.ap[ia][1],1.0])
        ia+=1
    newatoms.append(newatoms0)

for itype in range(poscar1.Ntype) :
    for newap in newatoms[itype] :
        poscar1.addatom(itype,newap)

for ia in range(poscar1.Natom) :
    poscar1.ap[ia][2]=(poscar1.ap[ia][2]-0.5)*poscar1.lc[2][2]/a3_new+0.5

poscar1.lc[2][2]=a3_new

poscar1.filewrite_vasp()
