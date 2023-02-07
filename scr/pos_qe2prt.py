#!/bin/python3

from dftscr.vaspfiles import poscar
import sys

if len(sys.argv)>1 :
    filename=sys.argv[1]
else :
    filename="scf.in"
poscar1=poscar.POSCAR(empty=True)
poscar1.fileread_qe(filename)
poscar1.filewrite_prt("prt.st")
