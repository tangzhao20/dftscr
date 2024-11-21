#!/usr/bin/env python3

# posmove.py a1 a2 a3
# move the vasp structure (Direct) by a1 a2 a3

# posmove.py a1 a2 a3 a4 a5 a6
# move the vasp structure (Direct) by a1/a2 a3/a4 a5/a6

from classes import Poscar
import sys

if len(sys.argv) < 4:
    print("too few argvs")
    print("shift by a1 a2 a3 or a1/a2 a3/a4 a5/a6")
    sys.exit()
elif len(sys.argv) < 7:
    disp = [float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])]
else:
    disp = [float(sys.argv[1])/float(sys.argv[2]), float(sys.argv[3]) /
            float(sys.argv[4]), float(sys.argv[5])/float(sys.argv[6])]

poscar1 = Poscar()
poscar1.read_vasp()

poscar1.move(disp=disp, lcart=False)

poscar1.write_vasp()
