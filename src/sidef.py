#!/usr/bin/env python3

# python3 sidef.py filename

# This code prints information of atoms in a silicon cluster
# A future version will allow creating the new structure with defect

import sys
import numpy as np
from classes import POSCAR

filename=sys.argv[1]
poscar1=POSCAR()
poscar1.fileread_vasp(filename)

if poscar1.atomtype[0]!="Si" :
    print("Error: The first atom is not Si.")
Nsi=poscar1.Naint[0]

poscar1.ap=(np.array(poscar1.ap)-0.5)
apc=poscar1.cartesian()[0:Nsi]
for i in range(Nsi) :
    d=(apc[i][0]**2+apc[i][1]**2+apc[i][2]**2)**0.5
    print(apc[i],d)

print(poscar1)
