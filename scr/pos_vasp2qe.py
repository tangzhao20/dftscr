#!/bin/python3

# Convert a POSCAR to QE structure format.

from dftscr.vaspfiles import poscar
import math

#f1=open("POSCAR","r")
#line=f1.readlines()
#f1.close()
#
#lc=[[],[],[]]
#lf=float(line[1].split()[0])
#for i in range(3) :
#    word0=line[i+2].split()
#    for j in range(3) :
#        lc[i].append(float(word0[j])*lf)
#satom=line[5].split()
#word0=line[6].split()
#Natom=0
#iatom=[]
#atom=[]
#for i in range(len(word0)) :
#    iatom.append(int(word0[i]))
#    Natom=Natom+iatom[i]
#for i in range(Natom) :
#    atom1=[]
#    word0=line[8+i].split()
#    for j in range(3) :
#        atom1.append(float(word0[j]))
#        if atom1[j]==-0.0 :
#            atom1[j]=0.0
#    atom.append(atom1)

poscar1=poscar.POSCAR()

#kgrid=[]
#for i in range(3) :
#    kgrid.append(math.floor(30.0/(poscar1.lc[i][0]**2+poscar1.lc[i][1]**2+poscar1.lc[i][2]**2)**0.5)+1)
#
#f2=open("qe.st","w")
#f2.write("CELL_PARAMETERS angstrom\n")
#for i in range(3) :
#    f2.write(f"  {poscar1.lc[i][0]:.12f}  {poscar1.lc[i][1]:.12f}  {poscar1.lc[i][2]:.12f}\n")
#f2.write("ATOMIC_SPECIES\n")
#for i in range(poscar1.Ntype) :
#    f2.write("  "+poscar1.atomtype[i]+"  1.0  PP\n")
#f2.write("ATOMIC_POSITIONS crystal\n")
#ij=0
#ik=0
#for i in range(poscar1.Natom) :
#    f2.write(f"  {poscar1.atomtype[ij]:2s}  {poscar1.ap[i][0]:.16f}  {poscar1.ap[i][1]:.16f}  {poscar1.ap[i][2]:.16f}\n")
#    ik=ik+1
#    if ik==poscar1.Naint[ij] :
#        ij=ij+1
#        ik=0
#f2.write("K_POINTS automatic\n")
#f2.write(f"  {kgrid[0]:d}  {kgrid[1]:d}  {kgrid[2]:d} 0 0 0\n")
#f2.close()

poscar1.filewrite_qe("qe.st")
