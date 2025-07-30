#!/usr/bin/env python3

# python3 posconvert.py package1 package2 (filename1)

# Convert a package1 to package2 structure format.

# Optional input: posconvert.in

import os
import sys
from classes import Poscar
from load_data import load_package_name

if len(sys.argv) < 3:
    print("Need input: package1 and package2")
    print("python3 posconvert.py package1 package2")
    sys.exit()

package_name = load_package_name()

package1 = sys.argv[1]
package2 = sys.argv[2]

# read filename1 from the command line
filename1 = ""
for iw in range(3, len(sys.argv)):
    if os.path.isfile(sys.argv[iw]):
        filename1 = sys.argv[iw]
        del sys.argv[iw]
        break

# read poscar
poscar1 = Poscar()
if package1 in package_name["vasp"]:
    poscar1.read_vasp(filename1)
elif package1 in package_name['qe']:
    poscar1.read_qe(filename1)
elif package1 in package_name['qexml']:
    poscar1.read_xml(filename1)
elif package1 in package_name['prt']:
    poscar1.read_prt(filename1)
elif package1 in package_name['parsec']:
    poscar1.read_parsec(filename1)
elif package1 in package_name['xyz']:
    poscar1.read_xyz(filename1)
else:
    print("Package "+package1+" input is not supported yet.")
    print("python3 posconvert.py package1 package2 (filename1)")
    sys.exit()

# read the Ndim
for iw in range(len(sys.argv)-1, -1, -1):
    if sys.argv[iw] in ["0d", "molecule", "cluster"]:
        poscar1.Ndim = 0
        del sys.argv[iw]
    elif sys.argv[iw] in ["1d", "wire"]:
        poscar1.Ndim = 1
        del sys.argv[iw]
    elif sys.argv[iw] in ["2d", "slab"]:
        poscar1.Ndim = 2
        del sys.argv[iw]
    elif sys.argv[iw] in ["3d", "bulk"]:
        poscar1.Ndim = 3
        del sys.argv[iw]

# read the file posconvert.in if it exists, then do some operations here
files = os.listdir()
if "posconvert.in" in files:
    f1 = open("posconvert.in", "r")
    line = f1.readlines()
    f1.close()
    for l in line:
        word = l.split()
        if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
            continue
        if word[0] == "wrap_to_cell":
            poscar1.wrap_to_cell()
        elif word[0] == "move":
            lcart = False
            if "cart" in word:
                lcart = True
                word.remove("cart")
            disp = [float(word[1]), float(word[2]), float(word[3])]
            poscar1.move(disp=disp, lcart=lcart)
        elif word[0] == "rotate":
            theta = float(word[1])
            poscar1.rotate(theta=theta)
        elif word[0] == "flip":
            poscar1.flip()
        elif word[0] == "supercell":
            N = [int(word[1]), int(word[2]), int(word[3])]
            poscar1.supercell(N=N)
        elif word[0] == "vacuum":
            z_vac = float(word[1])
            poscar1.vacuum(z_vac=z_vac)
        elif word[0] == "add_atom":
            new_type = word[1]
            ap = [float(word[2]), float(word[3]), float(word[4])]
            if new_type in poscar1.atomtype:
                itype = poscar1.atomtype.index(new_type)
                new_type = None
            else:
                itype = poscar1.Ntype
            poscar1.add_atom(itype=itype, ap=ap, new_type=new_type)
        elif word[0] == "delete_atom":
            ia = int(word[1]) - 1
            poscar1.delete_atom(ia=ia)
        elif word[0] == "replace_atom":
            ia = int(word[1]) - 1
            new_type = word[2]
            lhead = True
            if len(word) >= 4 and word[2].lower() == "tail":
                lhead = False
            poscar1.replace_atom(ia=ia, new_type=new_type, add_to_head=lhead)
        else:
            print("Warning: in posconvert.in, keyword "+word[0]+" is not supported yet")

# write poscar
if package2 in package_name["vasp"]:
    poscar1.write_vasp()
elif package2 in package_name['qe']:
    poscar1.write_qe()
elif package2 in package_name['prt']:
    poscar1.write_prt()
elif package2 in package_name['parsec']:
    lbohr = False
    lcart = False
    Ndim = -1
    for iw in range(len(sys.argv)-1, -1, -1):
        if sys.argv[iw].startswith("bohr") or sys.argv[iw].startswith("Bohr"):
            lbohr = True
            del sys.argv[iw]
        elif sys.argv[iw].startswith("cart") or sys.argv[iw].startswith("Cart"):
            lcart = True
            del sys.argv[iw]
    poscar1.write_parsec(lcartesian=lcart, lbohr=lbohr)
elif package2 in package_name['wannier90']:
    poscar1.write_wannier90()
elif package2 in package_name['xyz']:
    poscar1.write_xyz()
else:
    print("Package "+package2+" output is not supported yet.")
    print("python3 posconvert.py package1 package2")
    sys.exit()
