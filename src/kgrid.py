#!/usr/bin/env python3

"""
Generate k-grid in KPOINTS from POSCAR.

Usage: python3 kgrid.py
Input: POSCAR
Output: KPOINTS
"""

from classes import Poscar

poscar0 = Poscar()
poscar0.read_vasp()

k_grid = poscar0.find_k_grid()

f1 = open("KPOINTS", "w")
f1.write("Automatic mesh\n")
f1.write("0\n")
f1.write("Gamma\n")
f1.write(f"{k_grid[0]:d} {k_grid[1]:d} {k_grid[2]:d}\n")
f1.write("0 0 0\n")
f1.close()
