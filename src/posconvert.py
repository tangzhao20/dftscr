#!/usr/bin/env python3

# python3 posconvert.py package1 package2 (filename1)

# Convert a package1 to package2 structure format.

from classes import poscar
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
    if len(sys.argv)>=4 :
        filename1=sys.argv[3]
    else :
        filename1="POSCAR"
    poscar1=poscar.POSCAR(filename1)
elif package1 in packagename['qe'] :
    if len(sys.argv)>=4 :
        filename1=sys.argv[3]
    else :
        # find a scf.in or nscf.in file
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
    poscar1=poscar.POSCAR(empty=True)
    poscar1.fileread_qe(filename1)
elif package1 in packagename['prt'] :
    poscar1=poscar.POSCAR(empty=True)
    poscar1.fileread_prt("input")
elif package1 in packagename['qexml'] :
    poscar1=poscar.POSCAR(empty=True)
    poscar1.fileread_xml()
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
    poscar1.filewrite_parsec()
elif package2 in packagename['wannier90'] :
    poscar1.filewrite_wannier90()
else :
    print("Package "+package2+" output is not supported yet.")
    print("python3 posconvert.py package1 package2")
    sys.exit()

