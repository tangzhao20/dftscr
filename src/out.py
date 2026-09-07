#!/usr/bin/env python3

# Summarize the calculation results


import sys
from classes import Poscar, Outcar
from load_data import load_constant, load_package_name

if len(sys.argv) <= 1:
    print("python out.py package [prop1 prop2 ...]")
    sys.exit()
package = sys.argv[1]
custom_output = sys.argv[2:]

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
    print("Usage: python out.py package [prop1 prop2 ...]")
    sys.exit()

lengths, angles = poscar0.get_lattice_parameter()
volume = poscar0.volume()

toten_per_atom = outcar0.toten / poscar0.Natom

Js = outcar0.mag * mu_B / (volume * angstrom**3) * mu_0

if not custom_output:
    print(f"lattice parameters: {lengths[0]:.3f} {lengths[1]:.3f} {lengths[2]:.3f} A")
    print(f"lattice angles: {angles[0]:.1f} {angles[1]:.1f} {angles[2]:.1f} deg")
    print(f"volume: {volume:.3f} A^3")

    print()
    print(f"toten: {outcar0.toten:.6f} eV")
    print(f"toten per atom: {toten_per_atom:.6f} eV/atom")

    print()
    print(f"mag: {outcar0.mag:.6f} mu_B")
    print(f"J_s: {Js:.3f} T")

else:
    properties = {
        "a": lengths[0],
        "b": lengths[1],
        "c": lengths[2],
        "alpha": angles[0],
        "beta": angles[1],
        "gamma": angles[2],
        "volume": volume,
        "toten": outcar0.toten,
        "toten_per_atom": toten_per_atom,
        "mag": outcar0.mag,
        "js": Js,
    }

    aliases = {
        "vol": "volume",
        "e_per_atom": "toten_per_atom",
        "j_s": "js",
    }

    output_vals = []
    for item in custom_output:
        key = aliases.get(item.lower(), item.lower())
        if key in properties:
            val = properties[key]
            output_vals.append(f"{val:.6f}")
        else:
            output_vals.append("Unknown")

    print(" ".join(output_vals))
