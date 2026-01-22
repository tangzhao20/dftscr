#!/usr/bin/env python3

# Generate nscf k-point path for specific lattice
# python kgenerate.py (latt)
# Currently supported lattice: orcc, rhl

# input: POSCAR

import sys
import numpy as np
from classes import KpointsBand, Poscar
from load_data import load_constant

poscar0 = Poscar()
poscar0.read_vasp()

kpoints0 = KpointsBand()

if len(sys.argv) >= 2:
    latt = sys.argv[1].lower()
else:
    print("Error: please specify the lattice type")
    print("       currently supported: ORCC, RHL")
    sys.exit()

if latt.lower() == "orcc":
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

elif latt.lower() == "rhl":
    pi = load_constant("pi")
    degree = pi/180.0

    lengths = np.linalg.norm(poscar0.lc, axis=1)

    def get_angle(u, v):
        cos_theta = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        return np.arccos(cos_theta)

    angles = np.array([get_angle(poscar0.lc[1], poscar0.lc[2]),
                       get_angle(poscar0.lc[0], poscar0.lc[2]),
                       get_angle(poscar0.lc[0], poscar0.lc[1])])

    is_lengths_equal = np.allclose(lengths, lengths[0])
    is_angles_equal = np.allclose(angles, angles[0])

    if not (is_lengths_equal and is_angles_equal):
        print("Error: this is not a rhombohedral structure")
        print(f"  Lengths: {lengths} (A)")
        print(f"  Angles:  {angles/degree} deg")

    alpha = angles[0]
    a = lengths[0]
    if np.isclose(alpha, 0.5*pi):
        print("Error: alpha is 90 degrees. This is technically simple cubic, not RHL.")
        print(f"  a = {a:.6f} A, alpha = {alpha/degree:.6f} deg")
        sys.exit()

    if alpha < 0.5*pi:
        print(" Rhombohedral (RHL1, hR)")
        eta = (1.0 + 4 * np.cos(alpha)) / (2 + 4 * np.cos(alpha))
        nu = 0.75 - eta * 0.5

        kpoints0.kp_dict = {}
        kpoints0.kp_dict["Γ"] = [0.0, 0.0, 0.0]
        kpoints0.kp_dict["B"] = [eta, 0.5, 1.0-eta]
        kpoints0.kp_dict["B₁"] = [0.5, 1.0-eta, eta-1.0]
        kpoints0.kp_dict["F"] = [0.5, 0.5, 0.0]
        kpoints0.kp_dict["L"] = [0.5, 0.0, 0.0]
        kpoints0.kp_dict["L₁"] = [0.0, 0.0, -0.5]
        kpoints0.kp_dict["P"] = [eta, nu, nu]
        kpoints0.kp_dict["P₁"] = [1.0-nu, 1.0-nu, 1.0-eta]
        kpoints0.kp_dict["Q"] = [1.0-nu, nu, 0.0]
        kpoints0.kp_dict["X"] = [nu, 0.0, -nu]
        kpoints0.kp_dict["Z"] = [0.5, 0.5, 0.5]

        kpoints0.x_labels = [["Γ", "L", "B₁"], ["B", "Z", "Γ", "X"], ["Q", "F", "P₁", "Z"], ["L", "P"]]

    else:
        print(" Rhombohedral (RHL2, hR)")
        eta = 0.5 / np.tan(alpha * 0.5)
        nu = 0.75 - eta * 0.5

        kpoints0.kp_dict["Γ"] = [0.0, 0.0, 0.0]
        kpoints0.kp_dict["F"] = [0.5, -0.5, 0.0]
        kpoints0.kp_dict["L"] = [0.5, 0.0, 0.0]
        kpoints0.kp_dict["P"] = [1.0-nu, -nu, 1.0-nu]
        kpoints0.kp_dict["P₁"] = [nu, nu-1.0, nu-1.0]
        kpoints0.kp_dict["Q"] = [eta, eta, eta]
        kpoints0.kp_dict["Q₁"] = [1.0-eta, -eta, -eta]
        kpoints0.kp_dict["Z"] = [0.5, -0.5, 0.5]

        kpoints0.x_labels = [["Γ", "P", "Z", "Q", "Γ", "F", "P₁", "Q₁", "L", "Z"]]

    print(f"  a = {a:.6f} A, alpha = {alpha/degree:.6f} deg")

else:
    print("Error: lattice "+latt+" is not supported")
    print("       currently supported: ORCC, RHL")
    sys.exit()

kpoints0.write_kpathin()
