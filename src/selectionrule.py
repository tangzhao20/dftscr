#!/usr/bin/env python3

# This script calculates the selection rule of a point group, given the M matrix (character table of its irreps).

import sys
import os
import numpy as np
from fractions import Fraction

np.set_printoptions(formatter={'all': lambda x: str(Fraction(x).limit_denominator())})

point_group = sys.argv[1].lower()
file_name = point_group + ".dat"

this_dir, this_filename = os.path.split(__file__)
DATA_PATH = os.path.join(this_dir, "..", "data", "irreps", file_name)
f0 = open(DATA_PATH, "r")
line = f0.readlines()
f0.close()

M = []
irreps = []
for il in range(1, len(line)):
    word = line[il].split()
    if len(word) == 0 or word[0][0] in {"#", "!"}:
        continue
    M.append([])
    irreps.append(word[0])
    for iw in range(1, len(word)):
        M[-1].append(int(word[iw]))

M = np.array(M)
NM = len(M)
print(" NM = ", NM)
print(" irreps :")
print(irreps)
print(" M :")
print(M)
M_inv = np.linalg.inv(M)
# print("M_inv = ")
# print(M_inv)

# Here we print the nonzero integrals of <i|O|j>, where O is an isotropic operator
# Nonzero requires irreps[ii]==irreps[jj] in principle.
print(" Nonzero integrals of <i|O|j> :")
for ii in range(NM):
    ai = M[ii, :]
    for ij in range(ii, NM):
        aj = M[ij, :]
        aa = np.multiply(ai, aj) @ M_inv
        if abs(aa[0]) > 1.e-6:
            print(irreps[ii], irreps[ij])


# Here we print the nonzero integrals <ij|O|kl>, where O is an isotropic operator
print(" Nonzero integrals of <ij|O|kl> :")
count_all = 0
count_nonzero = 0
for ii in range(NM):
    ai = M[ii, :]
    for ij in range(ii, NM):
        aj = M[ij, :]
        aij = np.multiply(ai, aj)
        for ik in range(ij, NM):
            ak = M[ik, :]
            aijk = np.multiply(aij, ak)
            for il in range(ik, NM):
                al = M[il, :]
                aijkl = np.multiply(aijk, al)
                aa = aijkl @ M_inv
                count_all += 1
                if abs(aa[0]) > 1.e-6:
                    print(irreps[ii], irreps[ij], irreps[ik], irreps[il])
                    count_nonzero += 1
print(" NM^4 = ", NM**4)
print(" count_all = ", count_all)
print(" count_nonzero = ", count_nonzero)
