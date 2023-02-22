#!/usr/bin/env python3

# posmove.py a1 a2 a3
# move the vasp structure (Direct) by a1 a2 a3

# posmove.py a1 a2 a3 a4 a5 a6
# move the vasp structure (Direct) by a1/a2 a3/a4 a5/a6

from dftscr.vaspfiles import poscar
import sys

if len(sys.argv)<4 :
    print("too few argvs")
    print("shift by a1 a2 a3 or a1/a2 a3/a4 a5/a6")
    sys.exit()
elif len(sys.argv)<7 :
    factor=[float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3])]
else :
    factor=[float(sys.argv[1])/float(sys.argv[2]),float(sys.argv[3])/float(sys.argv[4]),float(sys.argv[5])/float(sys.argv[6])]

poscar1=poscar.POSCAR()
for i in range(poscar1.Natom) :
    for j in range(3) :
        poscar1.ap[i][j]=poscar1.ap[i][j]+factor[j]
poscar1.movetobox()

poscar1.filewrite("POSCAR.new")
