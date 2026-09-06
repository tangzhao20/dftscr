#!/usr/bin/env python3

# Summarize the calculation results


import sys
from classes import Poscar, Outcar
from load_data import load_constant, load_package_name

if len(sys.argv) <= 1:
    print("python out.py package")
    sys.exit()
package = sys.argv[1]

package_name = load_package_name()

mu_B = load_constant("mu_B")  # in J/T or A·m^2
mu_0 = load_constant("pi") * 4e-7
angstrom = load_constant("angstrom")

outcar0 = Outcar()
poscar0 = Poscar()

if package in package_name["vasp"]:
    # Inputs: CONTCAR, OUTCAR
    outcar0.read_vasp()
    poscar0.read_vasp("CONTCAR")
elif package in package_name["qe"]:
    # Input: *.xml
    outcar0.read_xml()
    poscar0.read_xml()
else:
    print("Package \""+package+"\" is not supported yet.")
    print("python out.py package")
    sys.exit()

lengths, angles = poscar0.get_lattice_parameter()
volume = poscar0.volume()

toten_per_atom = outcar0.toten / poscar0.Natom

Js = outcar0.mag * mu_B / (volume * angstrom**3) * mu_0

print(f"lattice parameters: {lengths[0]:.3f} {lengths[1]:.3f} {lengths[2]:.3f} A")
print(f"lattice angles: {angles[0]:.1f} {angles[1]:.1f} {angles[2]:.1f} deg")
print(f"volume: {volume:.3f} A^3")

print()
print(f"toten: {outcar0.toten:.6f} eV")
print(f"toten per atom: {toten_per_atom:.6f} eV/atom")

print()
print(f"mag: {outcar0.mag:.6f} mu_B")
print(f"J_s: {Js:.3f} T")

# print(f"{lengths[0]:.6f} {lengths[2]:.6f} {volume:.6f} {outcar0.toten:.6f} {outcar0.mag:.6f} {Js:.6f}")
