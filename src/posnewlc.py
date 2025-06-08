#!/usr/bin/env python3

# rotate the structure around the z-axis

# posnewlc.py

# input: posnewlc.in

from classes import Poscar
import sys

poscar1 = Poscar()
poscar1.read_vasp()

f0 = open("posnewlc.in", "r")
line = f0.readlines()
f0.close()

lc_new = []
for l in line:
    word = l.split()
    lc_new.append([float(word[0]), float(word[1]), float(word[2])])
poscar1.new_lc(lc_new)

poscar1.movetobox()

poscar1.write_vasp()
