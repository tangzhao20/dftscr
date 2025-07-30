#!/usr/bin/env python3

"""
Create structures with strains applied along x, y, or z direction.

Usage: python3 posstrain.py [x|y|z]
Input: POSCAR
Output: POSCAR_*
"""

import sys
import numpy as np
from classes import Poscar

N = 5
poscar0 = Poscar()
poscar0.read_vasp()

strain_range = [-0.1, 0.1]
strain_list = np.linspace(strain_range[0], strain_range[1], 2*N+1)

if len(sys.argv) < 2 or sys.argv[1] == "z":
    ix0 = 2
    str_dir = "z"
elif sys.argv[1] == "x":
    ix0 = 0
    str_dir = "x"
elif sys.argv[1] == "y":
    ix0 = 1
    str_dir = "y"
else:
    print("Error: Invalid argument. Please specify one of x, y, or z.\n")
    sys.exit(1)

print(f"Applying strain along the {str_dir} axis.")

lc_save = [0.0]*3
for ix1 in range(3):
    lc_save[ix1] = poscar0.lc[ix1, ix0]

for i in range(0, 2*N+1):
    for ix1 in range(3):
        poscar0.lc[ix1, ix0] = (1+strain_list[i]) * lc_save[ix1]
    poscar0.write_vasp("POSCAR_"+str(i))
