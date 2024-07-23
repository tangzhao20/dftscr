#!/usr/bin/env python3

# posflip.py
# flip the z coordinate of the structure

from classes import Poscar
import sys

poscar1=Poscar()
poscar1.read_vasp()

poscar1.flip()

poscar1.write_vasp()
