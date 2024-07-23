#!/usr/bin/env python3

# rotate the structure around the z-axis.

# posrotate.py theta(degrees)

from classes import Poscar
import sys

theta=float(sys.argv[1])

poscar1=Poscar()
poscar1.read_vasp()

poscar1.rotate(theta=theta)

poscar1.write_vasp()
