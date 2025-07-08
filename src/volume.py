#!/usr/bin/env python3

# volume.py
# print the volume of a structure

from classes import Poscar
import sys

poscar1 = Poscar()
poscar1.read_vasp()

vol = poscar1.volume()

print(f"volume: {vol:.6f} A^3")
