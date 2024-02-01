#!/usr/bin/env python3

# python3 posconvert.py package1 package2 (filename1)

# Convert a package1 to package2 structure format.

from classes import POSCAR
from commons import load_packagename
import sys
import os

if len(sys.argv)<3 :
    print("Need input: package1 and package2")
    print("python3 posconvert.py package1 package2")
    sys.exit()

packagename=load_packagename()

package1=sys.argv[1]
package2=sys.argv[2]

if package1 in packagename["vasp"] :
    filename1="POSCAR"
    for iw in range(3,len(sys.argv)) :
        if os.path.isfile(sys.argv[iw]) :
            filename1=sys.argv[iw]
            del sys.argv[iw]
            break
    poscar1=POSCAR(filename1)
elif package1 in packagename['qe'] :
    lname=False
    for iw in range(3,len(sys.argv)) :
        if os.path.isfile(sys.argv[iw]) :
            filename1=sys.argv[iw]
            lname=True
            del sys.argv[iw]
            break
    if lname==False :
        # find a scf/nscf/relax.in file
        files = os.listdir()
        if "scf.in" in files:
            filename1="scf.in"
        elif "nscf.in" in files:
            filename1="nscf.in"
        elif "relax.in" in files:
            filename1="relax.in"
        else :
            print("package1 is qe. Use input scf/nscf/relax.in, or add the filename at the end:")
            print("python3 posconvert.py qe package2 filename1")
            sys.exit()
    poscar1=POSCAR(empty=True)
    poscar1.fileread_qe(filename1)
elif package1 in packagename['qexml'] :
    poscar1=POSCAR(empty=True)
    poscar1.fileread_xml()
elif package1 in packagename['prt'] :
    poscar1=POSCAR(empty=True)
    poscar1.fileread_prt("input")
elif package1 in packagename['parsec'] :
    poscar1=POSCAR(empty=True)
    filename1="parsec.in"
    for iw in range(3,len(sys.argv)) :
        if os.path.isfile(sys.argv[iw]) :
            filename1=sys.argv[iw]
            del sys.argv[iw]
            break
    poscar1.fileread_parsec(filename=filename1)
elif package1 in packagename['xyz'] :
    filename1=""
    for iw in range(3,len(sys.argv)) :
        if os.path.isfile(sys.argv[iw]) :
            filename1=sys.argv[iw]
            del sys.argv[iw]
            break
    poscar1=POSCAR(empty=True)
    poscar1.fileread_xyz(filename1)
else :
    print("Package "+package1+" input is not supported yet.")
    print("python3 posconvert.py package1 package2")
    sys.exit()

if package2 in packagename["vasp"] :
    poscar1.filewrite()
elif  package2 in packagename['qe'] :
    poscar1.filewrite_qe()
elif package2 in packagename['prt'] :
    poscar1.filewrite_prt()
elif package2 in packagename['parsec'] :
    lbohr=False
    lcart=False
    Ndim=3
    for iw in range(len(sys.argv)-1,-1,-1) :
        if sys.argv[iw].startswith("molecule") or sys.argv[iw].startswith("cluster") :
            Ndim=0
            del sys.argv[iw]
        elif sys.argv[iw].startswith("slab") :
            Ndim=2
            del sys.argv[iw]
        elif sys.argv[iw].startswith("bohr") or sys.argv[iw].startswith("Bohr") :
            lbohr=True
            del sys.argv[iw]
        elif sys.argv[iw].startswith("cart") or sys.argv[iw].startswith("Cart") :
            lcart=True
            del sys.argv[iw]
    poscar1.filewrite_parsec(lcartesian=lcart,lbohr=lbohr,Ndim=Ndim)
elif package2 in packagename['wannier90'] :
    poscar1.filewrite_wannier90()
elif package2 in packagename['xyz'] :
    poscar1.filewrite_xyz()
else :
    print("Package "+package2+" output is not supported yet.")
    print("python3 posconvert.py package1 package2")
    sys.exit()

