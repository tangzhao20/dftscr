#!/usr/bin/env python3

# This script creates the structure of an h-BN monolayer (or flake) with or without a defect.

# Generate 2d h-BN structures by:
#  python3 hbndef.py N (defect)

# or generate h-BN flake structures by:
#  python3 hbndef.py flake r (defect)

# N defines the size of the supercell, r defines the radius of the h-BN flake.
# (defect) is the name of the defect in h-BN, such as CBVN.
# Use correct capitalization of atomic symbols.

import sys
import os
import math
import re
import numpy as np
from classes import Poscar

this_dir, this_filename = os.path.split(__file__)
hbn_structure = os.path.join(this_dir, "..", "data", "structures", "hbn.vasp")
poscar0 = Poscar()
poscar0.read_vasp(hbn_structure)
poscar0.Ndim = 2
bond_b_n = poscar0.lc[1][1]*2.0/3.0
bond_b_h = 1.20
bond_n_h = 1.02
# TODO: Here B-H bond and N-H bond are taken from wikipedia, with no reference.

if len(sys.argv) < 2:
    print("Need input:")
    print(" python3 hbndef.py N (defect)")
    print(" python3 hbndef.py flake r (defect)")
    sys.exit()

if "flake" in sys.argv:
    lflake = True
    sys.argv.remove("flake")
    for word in sys.argv:
        try:
            float(word)
            r_max = float(word)
            sys.argv.remove(word)
            break
        except ValueError:
            pass
    N = math.ceil(r_max / poscar0.lc[1][1]) * 2
    print("supercell: "+str(N)+", "+str(N)+", 1")
    poscar0.supercell([N, N, 1])
else:
    lflake = False
    for word in sys.argv:
        try:
            int(word)
            N = int(word)
            sys.argv.remove(word)
            break
        except ValueError:
            pass
    print("supercell: "+str(N)+", "+str(N)+", 1")
    poscar0.supercell([N, N, 1])
    poscar0.move([0.5, 0.5, 0])
ldef = False
if len(sys.argv) >= 2:  # contains a defect
    ldef = True

center = np.array([0.5, 0.5, 0.0]) @ np.array(poscar0.lc)  # center for Ndim = 2
apc = np.array(poscar0.cartesian())
apc_distance = np.linalg.norm(apc-center, axis=1)

if lflake:
    apc_b = apc[:N**2]
    apc_b_distance = apc_distance[:N**2]
    apc_n = apc[N**2:]
    apc_n_distance = apc_distance[N**2:]

    const = 3.0**0.5*0.5
    neighbor_b = np.array([[-const, -0.5, 0.0], [const, -0.5, 0.0], [0.0, 1.0, 0.0]])
    neighbor_n = np.array([[-const, 0.5, 0.0], [const, 0.5, 0.0], [0.0, -1.0, 0.0]])

    find_n = np.repeat(apc_b, 3, axis=0).reshape(-1, 3, 3) + neighbor_b*bond_b_n
    find_b = np.repeat(apc_n, 3, axis=0).reshape(-1, 3, 3) + neighbor_n*bond_b_n

    find_n_distance = np.linalg.norm(find_n-np.vstack([center] * 3), axis=2)
    find_b_distance = np.linalg.norm(find_b-np.vstack([center] * 3), axis=2)

    apc_h = []
    for ia0 in range(N**2-1, -1, -1):
        # N atoms
        ia = ia0 + N**2
        if apc_n_distance[ia0] > r_max:
            poscar0.delete_atom(ia)
        else:
            in_range = [0, 0, 0]
            for ia1 in range(3):
                if find_b_distance[ia0][ia1] < r_max + 1e-6:
                    in_range[ia1] = 1
            if sum(in_range) == 2:
                apc_h.append(apc_n[ia0] + neighbor_n[in_range.index(0)]*bond_n_h)
            if sum(in_range) == 1:
                poscar0.delete_atom(ia)
                apc_h.append(apc_n[ia0] + neighbor_n[in_range.index(1)]*(bond_b_n-bond_b_h))
    for ia0 in range(N**2-1, -1, -1):
        # B atoms
        if apc_b_distance[ia0] > r_max:
            poscar0.delete_atom(ia0)
        else:
            in_range = [0, 0, 0]
            for ia1 in range(3):
                if find_n_distance[ia0][ia1] < r_max + 1e-6:
                    in_range[ia1] = 1
            if sum(in_range) == 2:
                apc_h.append(apc_b[ia0] + neighbor_b[in_range.index(0)]*bond_b_h)
            if sum(in_range) == 1:
                poscar0.delete_atom(ia0)
                apc_h.append(apc_b[ia0] + neighbor_b[in_range.index(1)]*(bond_b_n-bond_n_h))

    ap_h = ((np.array(apc_h) @ np.linalg.inv(np.array(poscar0.lc))) + np.array([0.0, 0.0, 0.5])).tolist()
    new_type = "H"
    for ia in range(len(apc_h)):
        poscar0.add_atom(2, ap_h[ia], new_type=new_type, add_to_head=False)
        new_type = None

    apc = np.array(poscar0.cartesian())
    apc_distance = np.linalg.norm(apc-center, axis=1)

    radius = r_max + 5
    apc = apc - center + np.array([radius, radius, radius])
    center = np.array([radius, radius, radius])
    poscar0.ap = (apc*(1/radius/2.0)).tolist()
    poscar0.lc = [[radius*2.0, 0.0, 0.0], [0.0, radius*2.0, 0.0], [0.0, 0.0, radius*2.0]]
    poscar0.Ndim = 0

if ldef:  # create a defect

    defect = re.findall(r'[A-Z][a-z]?', sys.argv[1])
    if len(defect) not in [2, 4]:
        # Maybe this is not necessary?
        print("Error: Only single or double site defects are recognized here")
        sys.exit()
    defect_atoms = {}
    for i in range(len(defect)//2):
        defect_atoms[defect[i*2+1]] = defect[i*2]
    print(defect_atoms)

    b = np.argmin(apc_distance[:poscar0.Naint[0]])
    n = np.argmin(apc_distance[poscar0.Naint[0]:sum(poscar0.Naint[:2])])+poscar0.Naint[0]

    if "B" in defect_atoms:
        if defect_atoms["B"] == "V":
            poscar0.delete_atom(b)
            n -= 1
        else:
            poscar0.replace_atom(b, defect_atoms["B"], add_to_head=True)

    if "N" in defect_atoms:
        if defect_atoms["N"] == "V":
            poscar0.delete_atom(n)
        else:
            poscar0.replace_atom(n, defect_atoms["N"], add_to_head=True)

if not lflake:
    poscar0.move([-0.5, -0.5, 0])

filename = "hbn_"
if lflake:
    filename += "flake_" + str(r_max)
else:
    filename += str(N)
if ldef:
    filename += "_" + sys.argv[1]
filename += ".vasp"

poscar0.write_vasp(filename)

print(str(poscar0))
