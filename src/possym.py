#!/usr/bin/env python3

# Find the symmetry operations of a structure

import sys
import os
from classes import Poscar
from load_data import load_symops
from v3math import v3matchpp, v3dpp

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
    anew = []
    sumt = 0
    for it in range(poscar1.Ntype):
        anew0 = [-1]*poscar1.Naint[it]
        for ia1 in range(poscar1.Naint[it]):
            apc1 = apc[sumt+ia1]
            for ia2 in range(poscar1.Naint[it]):
                apc2 = apc[sumt+ia2]
                apc2_0 = [0.0]*3
                for ix1 in range(3):
                    for ix2 in range(3):
                        apc2_0[ix2] += mtx[isym][ix1][ix2]*apc2[ix1]
                if v3matchpp(apc1, apc2_0):
                    anew0[ia1] = ia2+sumt
                    break
        sumt += poscar1.Naint[it]
        anew.append(anew0)
    lwrite = True
    for anew_t in anew:
        for anew_a in anew_t:
            if anew_a == -1:
                lwrite = False
                break
        if lwrite == False:
            break
    if lwrite:
        # print(isym+1,ops[isym],mtx[isym])
        symlist.append(isym)
print(symlist)
for isym in symlist:
    print(str(isym)+" "+ops[isym]+" " + str(mtx[isym]))
