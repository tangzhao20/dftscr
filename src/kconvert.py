#!/usr/bin/env python3

# Convert a package1 to package2 k-path format
# python kconvert.py package1 package2 (N)

import sys
import os
from classes import KpointsBand, Poscar
from load_data import load_package_name

package_name = load_package_name()

Nk = 0  # Total Nk
Nk_line = 0  # same Nk on each line
if len(sys.argv) < 3:
    print("Need more arguments.")
    print("python kconvert.py package1 package2 (N)")
    sys.exit()
elif len(sys.argv) == 3:
    Nk_line = 40
else:
    for w in sys.argv[3:]:
        if w[0:3] == "Nk=":
            print(w)
            Nk = int(w[3:])
        elif w[0:8] == "Nk_line=":
            Nk_line = int(w[8:])
if Nk != 0 and Nk_line != 0:
    print("Only 1 of Nk or Nk_line should be set")
    sys.exit()

package1 = sys.argv[1]
package2 = sys.argv[2]

kpoints1 = KpointsBand()
if package1 in package_name["vasp"]:
    kpoints1.read_vasp()
elif package1 in package_name["kpathin"]:
    kpoints1.read_kpathin(Nk_line=Nk_line)
else:
    print("Package "+package1+" input is not supported yet.")
    print("python3 kconvert.py package1 package2 (Nk or Nk_line)")
    sys.exit()

if package2 in package_name["kpathin"]:
    kpoints1.write_kpathin()
elif package2 in package_name["qe"]:
    if Nk != 0:
        poscar1 = Poscar()
        poscar1.read_qe()
        rlc = poscar1.rlc()
        kpoints1.write_qe(Nk=Nk, rlc=rlc)
    else:
        kpoints1.write_qe()
elif package2 in package_name["vasp"]:
    kpoints1.write_vasp()
elif package2 in package_name["wannier90"]:
    kpoints1.write_wannier90()
elif package2 in package_name["parsec"]:
    kpoints1.write_parsec()
else:
    print("Package "+package2+" input is not supported yet.")
    print("python3 kconvert.py package1 package2 (N)")
    sys.exit()
