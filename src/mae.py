#!/usr/bin/env python3

import sys
import os
import time
import numpy as np
from classes import Poscar, Procar
from load_data import load_constant

start_time = time.time()

poscar0 = Poscar()
poscar0.read_vasp(filename="CONTCAR")

procar0 = Procar()
procar0.read_vasp()

Nb = procar0.Nb
for w in sys.argv:
    if w[0:3] == "Nb=":
        Nb = int(w[3:])
        print("Nb is reset to "+str(Nb))

# Identify the highest occupied state and the lowest hole state
occ_tol = 1e-3
occ_min = np.min(procar0.occ, axis=1)
occ_max = np.max(procar0.occ, axis=1)
cb = [0, 0]  # non-full bands: [cb:Nb]
vb = [0, 0]  # non-empty bands: [0:vb]
for ispin in range(2):
    cb[ispin] = np.where(occ_min[ispin] < 1.0 - occ_tol)[0][0]
    vb[ispin] = np.where(occ_max[ispin] > occ_tol)[0][-1] + 1
if Nb < max(vb):
    print("Warning: Nb is smaller than the highest occupied band.")

print(f"Read time: {time.time()-start_time:0.3f} s")
start_time = time.time()

atom_list = poscar0.atom_list()

# SOC factors taken from [New J. Phys. 21, 73054 (2019)]
soc_factor_dict = {"Fe": 0.05965, "Co": 0.07412}

atom_mask = np.isin(atom_list, list(soc_factor_dict))
soc_factors = np.array([soc_factor_dict[a] for a in np.array(atom_list)[atom_mask]])

# Lp = np.zeros((3, 3, 3))  # i: x, y, z
# # Order of orbs: py, pz, px
#
# Lp[0, 0, 1] = 1
# Lp[0, 1, 0] = -1
#
# Lp[1, 1, 2] = -1
# Lp[1, 2, 1] = 1
#
# Lp[1, 0, 2] = -1
# Lp[1, 2, 0] = 1

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
            e_diff = procar0.eig[ispin1, ik, 0:vb[ispin1]][:, np.newaxis] - \
                procar0.eig[ispin2, ik, cb[ispin2]:Nb][np.newaxis, :]
            e = e_diff / (e_diff**2 + eta**2)
            f = procar0.occ[ispin1, ik, 0:vb[ispin1]][:, np.newaxis] * \
                (1 - procar0.occ[ispin2, ik, cb[ispin2]:Nb][np.newaxis, :])
            # f = procar0.occ[ispin1, ik, :][:, np.newaxis] - procar0.occ[ispin2, ik, :][np.newaxis, :]
            e = e[np.newaxis, :, :]
            f = f[np.newaxis, :, :]

            c1 = procar0.complex[ispin1, ik, 0:vb[ispin1], atom_mask, 4:9].conj()  # numpy move the masked axis to front
            c2 = procar0.complex[ispin2, ik, cb[ispin2]:Nb, atom_mask, 4:9]

            # a: atoms; i,j: bands; x: directions; m,n: orbitals
            L = np.einsum("a, aim, xmn, ajn -> xij", soc_factors, c1, Ld, c2)
            e_i += sign * weight * np.sum(f * e * np.abs(L)**2 * 0.25, axis=(1, 2))

print("            x        y        z")
print(f"E_2PT   {e_i[0]*1e3:8.4f} {e_i[1]*1e3:8.4f} {e_i[2]*1e3:8.4f} meV")
e_i = e_i - np.min(e_i)
print(f"E_diff  {e_i[0]*1e3:8.4f} {e_i[1]*1e3:8.4f} {e_i[2]*1e3:8.4f} meV")

i_sort = np.argsort(e_i)
x_str = ["x", "y", "z"]
print(f"Med: {x_str[i_sort[1]]}, Easy: {x_str[i_sort[0]]}")

volume = poscar0.volume()
electron = load_constant("electron")
angstrom = load_constant("angstrom")
k1_med = e_i[i_sort[1]] / volume * electron / angstrom**3 / 1e6
print(f"K1_med-easy: {k1_med:0.3f} MJ/m3")
k1_hard = e_i[i_sort[2]] / volume * electron / angstrom**3 / 1e6
print(f"K1_hard-easy: {k1_hard:0.3f} MJ/m3")

print(f"Compute time: {time.time()-start_time:0.3f} s")
