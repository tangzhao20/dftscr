#!/usr/bin/env python3

# Calculate the effective mass of a phonon mode
# and use it as the weight to average the frequency

# Input: qpoints.yaml, POSCAR_i, POSCAR_f

import yaml
import matplotlib.pyplot as plt
from classes import Poscar
from pos2 import displacement
from load_data import load_constant

with open("qpoints.yaml") as f0:
    q1 = yaml.full_load(f0)

poscar_i = Poscar()
poscar_i.read_vasp(filename="POSCAR_i")
poscar_f = Poscar()
poscar_f.read_vasp(filename="POSCAR_f")

disp = displacement(poscar_i, poscar_f)

# read phonopy output file
Nb = len(q1["phonon"][0]["band"])
Na = len(q1["phonon"][0]["band"][0]["eigenvector"])
ph_freq = []  # freq[Nb]
ph_eigvec = []  # eigvec[Nb][Na][3]
for ib in range(Nb):
    ph_freq.append(q1["phonon"][0]["band"][ib]["frequency"])
    eigvec0 = []
    for ia in range(Na):
        eigvec1 = []
        for ix in range(3):
            eigvec1.append(q1["phonon"][0]["band"][ib]["eigenvector"][ia][ix][0])
        eigvec0.append(eigvec1)
    ph_eigvec.append(eigvec0)

overlap = [0.0]*Nb

for ib in range(Nb):
    for ia in range(Na):
        for j in range(3):
            overlap[ib] += ph_eigvec[ib][ia][j]*disp[ia][j]
    overlap[ib] = abs(overlap[ib])

freq_out = 0.0
summ = 0.0
for ib in range(Nb):
    freq_out += ph_freq[ib]*overlap[ib]
    summ += overlap[ib]
    print("%.3f" % ph_freq[ib]+" %.3f" % overlap[ib])

freq_out = freq_out/summ

print("average frequency: "+str(freq_out)+" THz")
planck = load_constant("planck")
electron = load_constant("electron")
hartree = load_constant("rydberg")*2
freq_out = freq_out*1e12*planck/electron/hartree
print(" = "+str(freq_out)+" a.u.")
