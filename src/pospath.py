#!/usr/bin/env python3

# Create structures as the linear combination of two initials.
# Input: POSCAR_i and POSCAR_f
# Output: POSCAR.*

# python3 pospath.py N (out)
# here N is the number of structures between initial two.

from classes import Poscar
from pos2 import match, moveatoms
import sys
import copy

lout=False
if "out" in sys.argv:
    lout=True
    sys.argv.remove("out")
N=int(sys.argv[1])
poscar0=Poscar()
poscar0.read_vasp(filename="POSCAR_i")
poscar1=Poscar()
poscar1.read_vasp(filename="POSCAR_f")
if match(poscar0, poscar1):
    pass
else:
    print("POSCARs don't match each other.")
    sys.exit()

if lout:
    for i in range(0,3*N+1):
        Ptmp=copy.deepcopy(poscar0) 
        moveatoms(Ptmp,poscar1,float(i)/float(N)-1.0)
        Ptmp.write_vasp("POSCAR."+str(i+1))
        del Ptmp
else:
    for i in range(0,N+1):
        Ptmp=copy.deepcopy(poscar0) 
        moveatoms(Ptmp,poscar1,float(i)/float(N))
        Ptmp.write_vasp("POSCAR."+str(i+1))
        del Ptmp
