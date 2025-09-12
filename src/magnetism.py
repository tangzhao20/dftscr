#!/usr/bin/env python3

# Script to calculate magnetic properties

# Inputs: CONTCAR, OUTCAR

from classes import Poscar, Outcar
from load_data import load_constant

outcar0 = Outcar()
outcar0.read_vasp()

poscar0 = Poscar()
poscar0.read_vasp("CONTCAR")

mu_B = load_constant("mu_B")  # in J/T or A·m^2
mu_0 = load_constant("pi") * 4e-7
angstrom = load_constant("angstrom")

Js = outcar0.mag * mu_B / (poscar0.volume() * angstrom**3) * mu_0

print(f"Js = {Js:.3f} T")
