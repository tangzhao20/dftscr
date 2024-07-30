#!/usr/bin/env python3

# python3 hbndef.py 6 (def)

# This code create the structure of h-BN monolayer (or flake) with a defect

import sys
import os
import numpy as np
from classes import Poscar

this_dir, this_filename = os.path.split(__file__)
hbn_structure = os.path.join(this_dir, "..", "data", "structures","hbn.vasp")
poscar0=Poscar()
poscar0.read_vasp(hbn_structure)
poscar0.Ndim=2

bond=poscar0.lc[1][1]*2/3
print(bond)

poscar0.supercell([6,6,1])
apc=poscar0.cartesian()
y_max=-1e6
y_min=1e6
for i in range(poscar0.Natom):
    if apc[i][1]>y_max :
        boron=i
        y_max=apc[i][1]
    if apc[i][1]<y_min :
        nitrogen=i
        y_min=apc[i][1]

print(boron,nitrogen)

poscar0.write_vasp("hbndef.vasp")

print(str(poscar0))

