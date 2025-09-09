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

soc_factor_dict = {"Fe": 0.05965, "Co": 0.07412}

atom_mask = np.isin(atom_list, list(soc_factor_dict))
soc_factors = np.array([soc_factor_dict[a] for a in np.array(atom_list)[atom_mask]])


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

eta = 0.0001

e_i = np.zeros(3)
for ik in range(procar0.Nk):
    weight = procar0.weight[ik]
    for ispin1 in range(2):
        for ispin2 in range(2):
            if ispin1 == ispin2:
                sign = 1
            else:
                sign = -1

            # vectorization over bands
            e_diff = procar0.eig[ispin1, ik, :][:, np.newaxis] - procar0.eig[ispin2, ik, :][np.newaxis, :]
            e = e_diff / (e_diff**2 + eta**2)
            f = procar0.occ[ispin1, ik, :][:, np.newaxis] * (1 - procar0.occ[ispin2, ik, :][np.newaxis, :])
            # f = procar0.occ[ispin1, ik, :][:, np.newaxis] - procar0.occ[ispin2, ik, :][np.newaxis, :]
            e = e[np.newaxis, :, :]
            f = f[np.newaxis, :, :]

            c1 = procar0.complex[ispin1, ik, :, atom_mask, 4:9].conj()  # numpy move the masked axis to front
            c2 = procar0.complex[ispin2, ik, :, atom_mask, 4:9]

            # a: atoms; i,j: bands; x: directions; m,n: orbitals
            L = np.einsum("a, aim, xmn, ajn -> xij", soc_factors, c1, Ld, c2)
            e_i += sign * weight * np.sum(f * e * np.abs(L)**2 * 0.25, axis=(1, 2))

            # for ib1 in range(procar0.Nb):  # higher band
            #     e1 = procar0.eig[ispin1, ik, ib1]
            #     occ1 = procar0.occ[ispin1, ik, ib1]
            #     for ib2 in range(procar0.Nb):  # lower band
            #         e2 = procar0.eig[ispin2, ik, ib2]
            #         occ2 = procar0.occ[ispin2, ik, ib2]
            #         #f = occ1
            #         f = occ1*(1-occ2)
            #         e = (e1-e2) / ((e1-e2)**2 + eta**2)
            #         L = np.zeros(3)
            #         for ia in range(poscar0.Natom):
            #             if atom_list[ia] == "B":
            #                 continue
            #             soc_factor = soc_factor_dict[atom_list[ia]]
            #             for ix in range(3):  # loop over x, y, z
            #                 L[ix] += np.abs(soc_factor * procar0.complex[ispin1, ik, ib1, ia, 4:9].conj() @ \
            #                          Ld[ix, :, :] @ procar0.complex[ispin2, ik, ib2, ia, 4:9])
            #         e_i += sign*weight*f*e*L**2
    print("ik: "+str(ik)+", e_i: "+str(e_i))

e_i = e_i - np.min(e_i)
print("\ne_i: "+str(e_i))
