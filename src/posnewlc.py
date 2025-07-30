#!/usr/bin/env python3

"""
Set the lattice vectors based on an input file

Usage: python3 posnewlc.py
Input: posnewlc.in
Output: POSCAR.new
"""

import sys
import numpy as np
from classes import Poscar

poscar1 = Poscar()
poscar1.read_vasp()

f0 = open("posnewlc.in", "r")
line = f0.readlines()
f0.close()

lc_new = np.zeros((3, 3))
for il in range(3):
    word = line[il].split()
    for ix in range(3):
        lc_new[il, ix] = float(word[ix])
poscar1.new_lc(lc_new)

poscar1.wrap_to_cell()

poscar1.write_vasp()
