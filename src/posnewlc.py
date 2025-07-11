#!/usr/bin/env python3

"""
Set the lattice vectors based on an input file

Usage: python3 posnewlc.py
Input: posnewlc.in
Output: POSCAR.new
"""

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
