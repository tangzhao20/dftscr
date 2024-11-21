#!/usr/bin/env python3

# from the conventional cell of silicon
# 1. create a supercell, whose size is determined by raidus
# 2. select the atoms inside of the sphere
# 3. add hydrogen atoms

# python3 silicon.py r

# where r is the radius in angstrom.

import sys
import os
import math
import numpy as np
from classes import Poscar

this_dir, this_filename = os.path.split(__file__)
silicon_structure = os.path.join(this_dir, "..", "data", "structures", "silicon.vasp")
poscar0 = Poscar()
poscar0.read_vasp(silicon_structure)
bond_si_si = poscar0.lc[0][1]*3.0**0.5/2.0
bond_si_h = 1.46  # need a reference

if len(sys.argv) < 2:
    print("Need input:")
    print(" python3 silicon.py r")
r_max = float(sys.argv[1])
N = math.ceil(r_max / (poscar0.lc[0][1]*2.0/3.0**0.5)) * 2
print("supercell: "+str(N)+", "+str(N)+", "+str(N))
poscar0.supercell([N, N, N])

center = np.array([0.5, 0.5, 0.5]) @ np.array(poscar0.lc)  # center for Ndim = 3
apc = np.array(poscar0.cartesian())
apc_distance = np.linalg.norm(apc-center, axis=1)

apc_si1 = apc[:N**3]
apc_si1_distance = apc_distance[:N**3]
apc_si2 = apc[N**3:]
apc_si2_distance = apc_distance[N**3:]

const = 3.0**(-0.5)
neighbor_si1 = np.array([[const, const, const], [const, -const, -const],
                        [-const, const, -const], [-const, -const, const]])
neighbor_si2 = np.array([[-const, -const, -const], [-const, const, const],
                        [const, -const, const], [const, const, -const]])

find_si2 = np.repeat(apc_si1, 4, axis=0).reshape(-1, 4, 3) + neighbor_si1*bond_si_si
find_si1 = np.repeat(apc_si2, 4, axis=0).reshape(-1, 4, 3) + neighbor_si2*bond_si_si

find_si2_distance = np.linalg.norm(find_si2-np.vstack([center] * 4), axis=2)
find_si1_distance = np.linalg.norm(find_si1-np.vstack([center] * 4), axis=2)

apc_h = []
for ia0 in range(N**3-1, -1, -1):
    # Si 2
    ia = ia0 + N**3
    if apc_si2_distance[ia0] > r_max:
        poscar0.delete_atom(ia)
    else:
        in_range = [0, 0, 0, 0]
        for ia1 in range(4):
            if find_si1_distance[ia0][ia1] < r_max + 1e-6:
                in_range[ia1] = 1
        if sum(in_range) == 2 or sum(in_range) == 3:
            for ia1 in range(4):
                if in_range[ia1] == 0:
                    apc_h.append(apc_si2[ia0] + neighbor_si2[ia1]*bond_si_h)
        if sum(in_range) == 1:
            poscar0.delete_atom(ia)
            apc_h.append(apc_si2[ia0] + neighbor_si2[in_range.index(1)]*(bond_si_si-bond_si_h))
for ia0 in range(N**3-1, -1, -1):
    # Si 1
    if apc_si1_distance[ia0] > r_max:
        poscar0.delete_atom(ia0)
    else:
        in_range = [0, 0, 0, 0]
        for ia1 in range(4):
            if find_si2_distance[ia0][ia1] < r_max + 1e-6:
                in_range[ia1] = 1
        if sum(in_range) == 2 or sum(in_range) == 3:
            for ia1 in range(4):
                if in_range[ia1] == 0:
                    apc_h.append(apc_si1[ia0] + neighbor_si1[ia1]*bond_si_h)
        if sum(in_range) == 1:
            poscar0.delete_atom(ia0)
            apc_h.append(apc_si1[ia0] + neighbor_si1[in_range.index(1)]*(bond_si_si-bond_si_h))

ap_h = (np.array(apc_h) @ np.linalg.inv(np.array(poscar0.lc))).tolist()
new_type = "H"
for ia in range(len(apc_h)):
    poscar0.add_atom(1, ap_h[ia], new_type=new_type, add_to_head=False)
    new_type = None

apc = np.array(poscar0.cartesian())
apc_distance = np.linalg.norm(apc-center, axis=1)

radius = r_max + 5
apc = apc - center + np.array([radius, radius, radius])
center = np.array([radius, radius, radius])
poscar0.ap = (apc*(1/radius/2.0)).tolist()
poscar0.lc = [[radius*2.0, 0.0, 0.0], [0.0, radius*2.0, 0.0], [0.0, 0.0, radius*2.0]]
poscar0.Ndim = 0

print(" radius "+str(r_max)+" A  Si "+str(poscar0.Naint[0])+" H "+str(poscar0.Naint[1]))
poscar0.write_vasp("Si"+str(poscar0.Naint[0])+"H"+str(poscar0.Naint[1])+".vasp")
