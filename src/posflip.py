#!/usr/bin/env python3

# posflip.py
# flip the z coordinate of the structure

from classes import POSCAR
import sys

poscar1=POSCAR()
poscar1.fileread_vasp()

for i in range(poscar1.Natom) :
    poscar1.ap[i][2]=1-poscar1.ap[i][2]
poscar1.movetobox()

poscar1.filewrite_vasp()
