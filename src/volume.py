#!/usr/bin/env python3

# volume.py (filename)
# print the volume of a structure

import sys
import os
from classes import Poscar

filename1 = "POSCAR"
for iw in range(1, len(sys.argv)):
    if os.path.isfile(sys.argv[iw]):
        filename1 = sys.argv[iw]
        # del sys.argv[iw]
        break

poscar1 = Poscar()
poscar1.read_vasp(filename=filename1)

vol = poscar1.volume()

print(f"volume: {vol:.12f} A^3")
