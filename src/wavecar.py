#!/usr/bin/env python3

# wavecar.py ik ib ispin
# output: wfn_ik_ib_ispin.vasp

import sys
from vaspwfc import vaspwfc
import os

ik = int(sys.argv[1])
ib = int(sys.argv[2])
ispin = int(sys.argv[3])

wav = vaspwfc('WAVECAR')
# try if the WAVECAR is calculated using Gamma version of VASP
try:
    phi = wav.wfc_r(ikpt=ik, iband=ib, ispin=ispin, ngrid=wav._ngrid * 2)
except:
    wav = vaspwfc('WAVECAR', lgamma=True)
    phi = wav.wfc_r(ikpt=ik, iband=ib, ispin=ispin, ngrid=wav._ngrid * 2)
wav.save2vesta(phi, poscar='POSCAR')

os.system("mv wfc_r.vasp wfn_"+str(ik)+"_"+str(ib)+"_"+str(ispin)+".vasp")
