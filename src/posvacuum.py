#!/usr/bin/env python3
#
# This script adds the vacuum layer to the structure.
#
# python posvacuum.py z_vac(A)
# where z_vac is in Angstrom.
#
# Input: POSCAR

from classes import Poscar
import sys

if len(sys.argv) < 2:
    print("Error: too few arguments")
    print("python posvacuum.py z_vac(A)")
    sys.exit()
else:
    z_vac = float(sys.argv[1])

poscar1 = Poscar()
poscar1.read_vasp()

poscar1.vacuum(z_vac)

poscar1.write_vasp()
