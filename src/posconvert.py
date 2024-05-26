#!/usr/bin/env python3

# python3 posconvert.py package1 package2 (filename1)

# Convert a package1 to package2 structure format.

# Optional input: posconvert.in

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

# read filename1 from the command line
filename1=""
for iw in range(3,len(sys.argv)) :
    if os.path.isfile(sys.argv[iw]) :
        filename1=sys.argv[iw]
        del sys.argv[iw]
        break

# read poscar
poscar1=POSCAR()
if package1 in packagename["vasp"] :
    if filename1=="" :
        filename1="POSCAR"
    poscar1.fileread_vasp(filename1)
elif package1 in packagename['qe'] :
    if filename1=="" :
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
    poscar1.fileread_qe(filename1)
elif package1 in packagename['qexml'] :
    poscar1.fileread_xml()
elif package1 in packagename['prt'] :
    poscar1.fileread_prt("input")
elif package1 in packagename['parsec'] :
    if filename1=="" :
        filename1="parsec.in"
    poscar1.fileread_parsec(filename=filename1)
elif package1 in packagename['xyz'] :
    poscar1.fileread_xyz(filename1)
else :
    print("Package "+package1+" input is not supported yet.")
    print("python3 posconvert.py package1 package2")
    sys.exit()

# read the file posconvert.in if it exists, then do some operations here
files = os.listdir()
if "posconvert.in" in files:
    f1=open("posconvert.in","r")
    line=f1.readlines()
    f1.close()
    for l in line :
        word=l.split()
        if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
            continue
        if word[0]=="movetobox" :
            poscar1.movetobox()
        elif word[0]=="move" :
            lcart=False
            if "cart" in word :
                lcart=True
                word.remove("cart")
            disp=[float(word[1]),float(word[2]),float(word[3])]
            poscar1.move(disp=disp,lcart=lcart)
        elif word[0]=="rotate" :
            theta=float(word[1])
            poscar1.rotate(theta=theta)
        elif word[0]=="flip" :
            poscar1.flip()
        elif word[0]=="supercell" :
            N=[int(word[1]),int(word[2]),int(word[3])]
            poscar1.supercell(N=N)
        elif word[0]=="vacuum" :
            z_vac=float(word[1])
            poscar1.vacuum(z_vac=z_vac)
        elif word[0]=="addatom" :
            itype=int(word[1])
            ap=[float(word[2]),float(word[3]),float(word[4])]
            newtype="X"
            if itype==poscar1.Ntype :
                newtype=word[5]
            poscar1.addatom(itype=itype,ap=ap,newtype=newtype)
        else :
            print("Warning: in posconvert.in, keyword "+word[0]+" is not supported yet")

# write poscar
if package2 in packagename["vasp"] :
    poscar1.filewrite_vasp()
elif  package2 in packagename['qe'] :
    poscar1.filewrite_qe()
elif package2 in packagename['prt'] :
    poscar1.filewrite_prt()
elif package2 in packagename['parsec'] :
    lbohr=False
    lcart=False
    Ndim=3
    for iw in range(len(sys.argv)-1,-1,-1) :
        if sys.argv[iw].startswith("molecule") or sys.argv[iw].startswith("cluster") or sys.argv[iw].startswith("0d") :
            Ndim=0
            del sys.argv[iw]
        elif sys.argv[iw].startswith("slab") or sys.argv[iw].startswith("2d") :
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

