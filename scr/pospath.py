#!/bin/python3
# Create structures as the linear combination of two initials.
# Input: POSCAR_i and POSCAR_f
# Output: POSCAR.*

# pospath.py N
# here N is the number of structures between initial two.

from dftscr.vaspfiles import poscar
import sys
import copy

N=int(sys.argv[1])
poscar0=poscar.POSCAR("POSCAR_i")
poscar1=poscar.POSCAR("POSCAR_f")
if poscar0.match(poscar1) :
    pass
else : 
    print("POSCARs don't match each other.")
    sys.exit()

for i in range(0,3*N+1) :
    Ptmp=copy.deepcopy(poscar0) 
    Ptmp.moveatoms(poscar1,float(i)/float(N)-1.0)
    Ptmp.filewrite("POSCAR."+str(i+1))
    del Ptmp
