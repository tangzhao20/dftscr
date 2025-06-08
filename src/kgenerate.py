#!/usr/bin/env python3

# Generate nscf k-point path for specific lattice
# python kgenerate.py (latt)
# Currently supported lattice: orcc

# input: POSCAR

import sys
import os
from classes import KpointsBand, Poscar
from load_data import load_package_name

poscar0 = Poscar()
poscar0.read_vasp()

kpoints0 = KpointsBand()

if len(sys.argv) >= 2:
    latt = sys.argv[1].lower()
else:
    print("Error: please specify the lattice type")
    print("       currently supported: orcc")
    sys.exit()

if latt == "orcc":
    print(" C-centered orthorhombic (ORCC, oS)")
    print("  a1 = (a/2, -b/2,  0)")
    print("  a2 = (a/2,  b/2,  0)")
    print("  a3 = (  0,    0,  c)")
    print("  where a < b")

    a = poscar0.lc[0][0] * 2.0
    b = - poscar0.lc[0][1] * 2.0
    if a >= b:
        print("Error: a < b is required")
        sys.exit()

    zeta = (1.0 + (a/b)**2) / 4.0

    kpoints0.kpdict = {}
    kpoints0.kpdict["Γ"] = [0.0, 0.0, 0.0]
    kpoints0.kpdict["A"] = [zeta, zeta, 0.5]
    kpoints0.kpdict["A₁"] = [-zeta, 1.0-zeta, 0.5]
    kpoints0.kpdict["R"] = [0.0, 0.5, 0.5]
    kpoints0.kpdict["S"] = [0.0, 0.5, 0.0]
    kpoints0.kpdict["T"] = [-0.5, 0.5, 0.5]
    kpoints0.kpdict["X"] = [zeta, zeta, 0.0]
    kpoints0.kpdict["X₁"] = [-zeta, 1.0-zeta, 0.0]
    kpoints0.kpdict["Y"] = [-1/2, 1/2, 0.0]
    kpoints0.kpdict["Z"] = [0.0, 0.0, 0.5]

    kpoints0.xlabels = [["Γ", "X", "S", "A", "Z", "Γ", "Y", "X₁", "A₁", "T", "Y"], ["Z", "T"]]

else:
    print("Error: lattice "+latt+" is not supported")
    print("       currently supported: orcc")
    sys.exit()

kpoints0.write_kpathin()
