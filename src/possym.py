#!/usr/bin/env python3

# Find the symmetry operations of a structure

import sys
import os
import numpy as np
from classes import Poscar
from load_data import load_symops


# The structure is read from parsec format
lname = "parsec.in"
for iw in range(1, len(sys.argv)):
    if os.path.isfile(sys.argv[iw]):
        filename1 = sys.argv[iw]
        del sys.argv[iw]
        break
poscar1 = Poscar()
poscar1.read_parsec(filename1)

mtx, ops = load_symops()

apc = poscar1.cartesian()

symlist = []
for isym in range(48):
    atom_match = []
    apc_new = apc @ mtx[isym, :, :].T
    sum_type = 0
    for it in range(poscar1.Ntype):
        atom_match0 = [-1]*poscar1.Naint[it]
        for ia1 in range(poscar1.Naint[it]):
            for ia2 in range(poscar1.Naint[it]):
                if np.allclose(apc[sum_type+ia1, :], apc_new[sum_type+ia2, :]):
                    atom_match0[ia1] = ia2+sum_type
                    break
        sum_type += poscar1.Naint[it]
        atom_match = atom_match + atom_match0
    lwrite = True
    if -1 not in atom_match:
        symlist.append(isym)

print(symlist)
for isym in symlist:
    print(f"{isym:2d} {ops[isym]:9s} ", end="")
    for ix1 in range(3):
        for ix2 in range(3):
            mtx_out = int(round(mtx[isym, ix1, ix2]))
            print(f"{mtx_out:2d} ", end="")
        print(" ", end="")
    print("")
