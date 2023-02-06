#!/bin/python3

# create nscf k-point path for QE band structure calculations

# input: kpath.in

import sys
from dftscr.vaspfiles import kpoints_band

kpoints1=kpoints_band.KPOINTS_band(empty=True)

if len(sys.argv)==1 :
    N=40
else :
    N=int(sys.argv[1])

kpoints1.fileread_kpathin(nk_line=N)

kpoints1.filewrite_qe()
