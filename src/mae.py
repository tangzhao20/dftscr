#!/usr/bin/env python3

import sys
import os
import numpy as np
from classes import Poscar, Procar

poscar0 = Poscar()
poscar0.read_vasp(filename="CONTCAR")

procar0 = Procar()
procar0.read_vasp()

atom_list = poscar0.atom_list()

eta = 0.0001

soc_factor_dict = {"Fe": 0.05965, "Co": 0.07412, "B": 0.0}

Lp = np.zeros((3, 3, 3))  # i: x, y, z
# Order of orbs: py, pz, px

Lp[0, 0, 1] = 1
Lp[0, 1, 0] = -1

Lp[1, 1, 2] = -1
Lp[1, 2, 1] = 1

Lp[1, 0, 2] = -1
Lp[1, 2, 0] = 1

Ld = np.zeros((3, 5, 5))  # i: x, y, z
# Order of orbs: dxy, dyz, dz2, dxz, x2-y2

Ld[0, 0, 3] = 1.0
Ld[0, 1, 2] = 3.0**0.5
Ld[0, 1, 4] = 1.0
Ld[0, 2, 1] = -3.0**0.5
Ld[0, 3, 0] = -1.0
Ld[0, 4, 1] = -1.0

Ld[1, 0, 1] = 1.0
Ld[1, 1, 0] = -1.0
Ld[1, 2, 3] = -3.0**0.5
Ld[1, 3, 2] = 3.0**0.5
Ld[1, 3, 4] = -1.0
Ld[1, 4, 3] = 1.0

Ld[2, 0, 4] = -2.0
Ld[2, 1, 3] = -1.0
Ld[2, 3, 1] = 1.0
Ld[2, 4, 0] = 2.0

e_i = np.zeros(3)
L = np.zeros(3)
for ik in range(procar0.Nk):
    weight = procar0.weight[ik]
    for ia in range(poscar0.Natom):
        soc_factor = soc_factor_dict[atom_list[ia]]
        print("ik: "+str(ik)+", weight: "+str(weight)+", ia: "+str(ia) +
              ", atom: "+str(atom_list[ia])+", soc_factor: "+str(soc_factor))
        for ispin1 in range(2):
            for ispin2 in range(2):
                if ispin1 == ispin2:
                    sign = 1
                else:
                    sign = -1
                for ib1 in range(procar0.Nb):  # higher band
                    e1 = procar0.eig[ispin1, ik, ib1]
                    occ1 = procar0.occ[ispin1, ik, ib1]
                    for ib2 in range(procar0.Nb):  # lower band
                        e2 = procar0.eig[ispin2, ik, ib2]
                        occ2 = procar0.occ[ispin2, ik, ib2]
                        # f = occ1 - occ2
                        f = occ1*(1-occ2)
                        e = (e1-e2) / ((e1-e2)**2 + eta**2)
                        for ix in range(3):  # loop over x, y, z
                            L[ix] = procar0.proj[ispin1, ik, ib1, ia, 4:9] @ Ld[ix,
                                                                                :, :] @ procar0.proj[ispin2, ik, ib2, ia, 4:9]
                            e_i[ix] += sign*weight*soc_factor**2*f*e*L[ix]**2
                        # print("ispin1: "+str(ispin1)+", ispin2: "+str(ispin2)+", sign: "+str(sign)+
                        #      ", ib1: "+str(ib1)+", ib2: "+str(ib2)+", f: "+str(f)+", e: "+str(e))
                        # print("L: "+str(L)+", e_i: "+str(e_i))
        print("e_i: "+str(e_i))
        print()

print("e_i: "+str(e_i))
