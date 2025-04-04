#!/usr/bin/env python3

# Construct the supercell from unit cell

# supercell.py n1 n2 n3

# input: POSCAR
# output: POSCAR.new

from classes import Poscar
import sys

poscar1 = Poscar()
poscar1.read_vasp()

if len(sys.argv) < 4:
    print("too few argv")
    print("expand by n1 n2 n3")
    sys.exit()
N = [int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])]

poscar1.supercell(N)

poscar1.write_vasp()
