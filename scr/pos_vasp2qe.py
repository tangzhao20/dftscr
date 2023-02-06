#!/bin/python3

# Convert a POSCAR to QE structure format.

from dftscr.vaspfiles import poscar
import math

poscar1=poscar.POSCAR()

poscar1.filewrite_qe("qe.st")
