#!/usr/bin/env python3

# create nscf k-point path for QE band structure calculations
# python kconvert.py package1 package2 (N)

# input: kpath.in

import sys
import os
from classes import KPOINTS_band, POSCAR
from load_data import load_packagename

packagename=load_packagename()

Nk=0 # Total Nk
Nk_line=0 # same Nk on each line
if len(sys.argv)<3:
    print("Need more arguments.")
    print("python kconvert.py package1 package2 (N)")
    sys.exit()
elif len(sys.argv)==3 :
    Nk_line=40
else :
    for w in sys.argv[3:] :
        if w[0:3]=="Nk=" :
            print(w)
            Nk=int(w[3:])
        elif w[0:8]=="Nk_line=" :
            Nk_line=int(w[8:])
if Nk!=0 and Nk_line!=0 :
    print("Only 1 of Nk or Nk_line should be set")
    sys.exit()
    
package1=sys.argv[1]
package2=sys.argv[2]

kpoints1=KPOINTS_band()
if package1 in packagename["vasp"] :
    kpoints1.fileread_vasp()
elif package1 in packagename["kpathin"] :
    kpoints1.fileread_kpathin(Nk_line=Nk_line)
else :
    print("Package "+package1+" input is not supported yet.")
    print("python3 kconvert.py package1 package2 (Nk or Nk_line)")
    sys.exit()

if package2 in packagename["kpathin"] :
    kpoints1.filewrite_kpathin()
elif package2 in packagename["qe"] :
    if Nk!=0 :
        poscar1=POSCAR()
        poscar1.fileread_qe()
        rlc=poscar1.rlc()
        kpoints1.filewrite_qe(Nk=Nk,rlc=rlc)
    else :
        kpoints1.filewrite_qe()
elif package2 in packagename["vasp"] :
    kpoints1.filewrite_vasp()
elif package2 in packagename["wannier90"] :
    kpoints1.filewrite_wannier90()
elif package2 in packagename["parsec"] :
    kpoints1.filewrite_parsec()
else :
    print("Package "+package2+" input is not supported yet.")
    print("python3 kconvert.py package1 package2 (N)")
    sys.exit()
