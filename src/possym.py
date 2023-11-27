#!/usr/bin/env python3

# Find the symmetry operations of a structure

import sys
import os
from classes import POSCAR

from commons import load_symops
from v3math import v3matchpp, v3dpp

# The structure is read from parsec format
lname=False
for iw in range(1,len(sys.argv)) :
    if os.path.isfile(sys.argv[iw]) :
        filename1=sys.argv[iw]
        lname=True
        del sys.argv[iw]
        break
poscar1=POSCAR(empty=True)
if lname :
    poscar1.fileread_parsec(filename1)
else :
    poscar1.fileread_parsec()

mtx, ops = load_symops()

# Here I use the direct coordinate, since this package assume cubic cell 
# for cluster structures in the format of parsec.

ap=[]
for ia in range(poscar1.Natom) :
    ap0=[0.0]*3
    for ix1 in range(3) :
        for ix2 in range(3) :
            ap0[ix2]+=poscar1.lc[ix1][ix2]*(poscar1.ap[ia][ix1]-0.5)
    ap.append(ap0)

symlist=[]
for isym in range(48) :
    anew=[]
    sumt=0
    for it in range(poscar1.Ntype) :
        anew0=[-1]*poscar1.Naint[it]
        for ia1 in range(poscar1.Naint[it]) :
            ap1=ap[sumt+ia1]
            for ia2 in range(poscar1.Naint[it]) :
                ap2=ap[sumt+ia2]
                ap2_0=[0.0]*3
                for ix1 in range(3) :
                    for ix2 in range(3) :
                        ap2_0[ix2]+=mtx[isym][ix1][ix2]*ap2[ix1]
                if v3matchpp(ap1, ap2_0) :
                    anew0[ia1]=ia2+sumt
                    break
        sumt+=poscar1.Naint[it]
        anew.append(anew0)
    lwrite=True
    for anew_t in anew :
        for anew_a in anew_t :
            if anew_a==-1 :
                lwrite=False
                break
        if lwrite==False :
            break
    if lwrite :
        #print(isym+1,ops[isym],mtx[isym])
        symlist.append(isym+1)
print(symlist)
