#!/usr/bin/env python3

# Move the structure according to a phonon mode. 
# Typically designed for a JTE style broken symmetry of a point defect by soft modes.

# phmove.py ib
# where ib is the index of phonon mode (start from 1)

import sys
import numpy as np
import yaml
from classes import Poscar

if len(sys.argv)==1:
    print("No index of phonon mode. set to 1.")
    ib=1
else:
    ib=int(sys.argv[1])

poscar0=Poscar()
poscar0.read_vasp()

with open("qpoints.yaml") as f0:
    q1=yaml.full_load(f0)

Na=q1["natom"]
if Na!=poscar0.Natom :
    print("Error: Number of atoms from POSCAR and qpoints.yaml doesn`t match")
    sys.exit()

Nb=len(q1["phonon"][0]["band"])
if ib>Nb :
    print("Error: Index of phonon mode > number of phonon modes.")
    sys.exit()

ph_eigvec=[]
for ia in range(Na) :
    eigvec0=[]
    for ix in range(3) :
        eigvec0.append(q1["phonon"][0]["band"][ib-1]["eigenvector"][ia][ix][0])
    ph_eigvec.append(eigvec0)

ph_eigvec=np.array(ph_eigvec)
disp=ph_eigvec/np.sum(np.linalg.norm(ph_eigvec,axis=1))*0.4 # scale: 0.4 A

poscar0.move(disp=disp,lcart=True)

poscar0.write_vasp("POSCAR."+str(ib))
