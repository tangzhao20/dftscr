#!/usr/bin/env python3
#
# This script calculates the distance between two atoms.
# Note: the periodicity is not considered yet and should be fixed in a future version.
#
# python distance.py ia1 ia2 filename
#
# Input: *.vasp

from classes import Poscar
import sys
import os
import numpy as np

filename1 = ""
for iw in range(1, len(sys.argv)):
    if os.path.isfile(sys.argv[iw]):
        filename1 = sys.argv[iw]
        del sys.argv[iw]
        break

ia1 = int(sys.argv[1])
ia2 = int(sys.argv[2])

poscar1 = Poscar()
poscar1.read_vasp(filename=filename1)

apc = np.array(poscar1.cartesian())

dis = np.linalg.norm(apc[ia1, :] - apc[ia2, :])

print(f"{dis:10.5f}")
