#!/usr/bin/env python3

# posflip.py
# flip the z coordinate of the structure

from classes import POSCAR
import sys

poscar1=POSCAR()
poscar1.fileread_vasp()

poscar1.flip()

poscar1.filewrite_vasp()
