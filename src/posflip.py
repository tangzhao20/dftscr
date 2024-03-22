#!/usr/bin/env python3

# posflip.py
# flip the z coordinate of the structure

from classes import POSCAR
import sys

poscar1=POSCAR()
for i in range(poscar1.Natom) :
    for j in range(3) :
        poscar1.ap[i][j]=1-poscar1.ap[i][j]
poscar1.movetobox()

poscar1.filewrite("POSCAR.new")
