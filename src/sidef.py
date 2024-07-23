#!/usr/bin/env python3

# python3 sidef.py filename

# This code prints information of atoms in a silicon cluster
# A future version will allow creating the new structure with defect

import sys
import numpy as np
from copy import deepcopy
from classes import Poscar

filename=sys.argv[1]
poscar1=Poscar()
poscar1.read_vasp(filename)
poscar1.Ndim=0

if poscar1.atomtype[0]!="Si" :
    print("Error: The first atom is not Si.")
Nsi=poscar1.Naint[0]

apc=poscar1.cartesian()[0:Nsi]

print(poscar1)

def info_defect(i) :
    d=(apc[i][0]**2+apc[i][1]**2+apc[i][2]**2)**0.5
    poscar_new=deepcopy(poscar1)
    poscar_new.delete_atom(i)
    output_name="parsec_st_si"+str(poscar1.Naint[0])+"_"+str(i+1)+".dat"
    poscar_new.write_parsec(filename=output_name)
    print(f"{i+1:5d} {apc[i][0]:9.4f} {apc[i][1]:9.4f} {apc[i][2]:9.4f} {d:9.4f}")

tol=1e-4

print("center")
for i in range(Nsi) :
    if abs(apc[i][0]) < tol and abs(apc[i][1]) < tol and abs(apc[i][2]) < tol :
        i_save=i
        info_defect(i)
print()

print("[111]")
for i in range(Nsi) :
    if abs(apc[i][0]-apc[i][1]) < tol and abs(apc[i][0]-apc[i][2]) < tol :
        if i==i_save :
            continue
        info_defect(i)
print()

print("[110]")
for i in range(Nsi) :
    if abs(apc[i][0]-apc[i][1]) < tol and abs(apc[i][2]) < tol :
        if i==i_save :
            continue
        info_defect(i)
print()

print("[100]")
for i in range(Nsi) :
    if abs(apc[i][1]) < tol and abs(apc[i][2]) < tol :
        if i==i_save :
            continue
        info_defect(i)
print()
