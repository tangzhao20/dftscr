#!/usr/bin/env python3

# Print the toten vs. distortion coordinate

# Input: POSCAR_i, pos_*/OSZICAR, pos_*/POSCAR, (pos_*/EIGENVAL)

from classes import Poscar, Eigenval
from pos2 import general_q

N = 41

# True: total energy
# False: eigenvalues
is_toten = True

poscar_0 = Poscar()
poscar_0.read_vasp(filename="POSCAR_i")

energy = []
pos = []

for i in range(N):

    if is_toten:
        f1 = open(str(i+1)+"/OSZICAR")
        line = f1.readlines()
        f1.close()
        word = line[-1].split()
        energy.append(float(word[2]))
    else:
        # not needed
        eigenval_1 = Eigenval()
        eigenval_1.read_vasp(str(i+1)+"/EIGENVAL")
        energy0 = []
        for ispin in range(2):
            for ib in range(140, 145):
                energy0.append(eigenval_1.eig[0][ib][ispin])
        energy.append(energy0)
        del energy0

    poscar_1 = Poscar()
    poscar_1.read_vasp(filename=str(i+1)+"/POSCAR")

# Q = Renormalization factor:
#    pos.append(poscar_1.total_distance(poscar_0))

#   Q = Generalized coordinates (renormalization factor weighted by atomic mass)
    pos.append(general_q(poscar_0, poscar_1))

# Q = Displacement of the central factor
#    pos.append(poscar_1.lc[2, 2]*poscar_1.ap[35][2])

    del poscar_1

if is_toten:
    energy_zero = energy[0]
    for i in range(N):
        energy[i] = energy[i]-energy_zero
else:
    energy_zero = energy[0][2]
    for i in range(N):
        for j in range(len(energy[i])):
            energy[i][j] = energy[i][j]-energy_zero

energy_out = []
pos_out = []
for i in range(N-1):
    pos_out.append(-pos[-i-1])
    energy_out.append(energy[-i-1])

for i in range(N):
    pos_out.append(pos[i])
    energy_out.append(energy[i])

if is_toten:
    fout = open("toten_fit.dat", "w")
    fout.write("Q toten(eV)\n")
    for i in range(2*N-1):
        fout.write("%.6f" % pos_out[i]+" %.6f\n" % energy_out[i])
    fout.close()
else:
    # not needed
    fout = open("eig_fit.dat", "w")
    for i in range(2*N-1):
        fout.write("%.6f" % pos_out[i])
        for j in range(len(energy_out[i])):
            fout.write(" %.6f" % energy_out[i][j])
        fout.write("\n")
    fout.close()
