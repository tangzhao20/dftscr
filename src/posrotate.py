#!/usr/bin/env python3

# rotate the structure around the z-axis.

# posrotate.py theta(degrees)

from classes import POSCAR
import sys

theta=float(sys.argv[1])

poscar1=POSCAR()
poscar1.fileread_vasp()

poscar1.rotate(theta=theta)

poscar1.filewrite_vasp()
