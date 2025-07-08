#!/usr/bin/env python3

# Create structures with strains on z direction
# Input: POSCAR
# Output: POSCAR_*

# python3 posstrain.py

import sys
import numpy as np
from classes import Poscar

N = 5
poscar0 = Poscar()
poscar0.read_vasp()

strain_range = [-0.1, 0.1]
strain_list = np.linspace(strain_range[0], strain_range[1], 2*N+1)
print(strain_list)

if max(abs(poscar0.lc[2][0]), abs(poscar0.lc[2][1])) > 1.e-6:
    print("Warning: lattice vector a2 is not vertical")
z0 = poscar0.lc[2][2]

for i in range(0, 2*N+1):
    poscar0.lc[2][2] = (1+strain_list[i]) * z0
    poscar0.write_vasp("POSCAR_"+str(i))
