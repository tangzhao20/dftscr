#!/usr/bin/env python3

# Generate nscf k-point path for specific lattice
# python kgenerate.py (latt)
# Currently supported lattice: orcc

# input: POSCAR

import sys
from classes import KpointsBand, Poscar

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

    a = poscar0.lc[0, 0] * 2.0
    b = - poscar0.lc[0, 1] * 2.0
    if a >= b:
        print("Error: a < b is required")
        sys.exit()

    zeta = (1.0 + (a/b)**2) / 4.0

    kpoints0.kp_dict = {}
    kpoints0.kp_dict["Γ"] = [0.0, 0.0, 0.0]
    kpoints0.kp_dict["A"] = [zeta, zeta, 0.5]
    kpoints0.kp_dict["A₁"] = [-zeta, 1.0-zeta, 0.5]
    kpoints0.kp_dict["R"] = [0.0, 0.5, 0.5]
    kpoints0.kp_dict["S"] = [0.0, 0.5, 0.0]
    kpoints0.kp_dict["T"] = [-0.5, 0.5, 0.5]
    kpoints0.kp_dict["X"] = [zeta, zeta, 0.0]
    kpoints0.kp_dict["X₁"] = [-zeta, 1.0-zeta, 0.0]
    kpoints0.kp_dict["Y"] = [-1/2, 1/2, 0.0]
    kpoints0.kp_dict["Z"] = [0.0, 0.0, 0.5]

    kpoints0.x_labels = [["Γ", "X", "S", "A", "Z", "Γ", "Y", "X₁", "A₁", "T", "Y"], ["Z", "T"]]

else:
    print("Error: lattice "+latt+" is not supported")
    print("       currently supported: orcc")
    sys.exit()

kpoints0.write_kpathin()
