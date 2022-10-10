#!/bin/python3

# convert a POSCAR to a paratec structure file

from dftscr.vaspfiles import poscar 

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
#
#volume=(lc[0][0]*(lc[1][1]*lc[2][2]-lc[1][2]*lc[2][1])+lc[0][1]*(lc[1][2]*lc[2][0]-lc[1][0]*lc[2][2])+lc[0][2]*(lc[1][0]*lc[2][1]-lc[1][1]*lc[2][0]))*6.748334503468005

poscar1=poscar.POSCAR()

#volume=poscar1.volume()*6.748334503468005
## convert Angstrom^3 to bohr^3
#
#f2=open("input.coord","w")
#f2.write("begin latticevecs\n")
#for i in range(3) :
#    f2.write("coord  "+f"{poscar1.lc[i][0]:20.16f}"+"  "+f"{poscar1.lc[i][1]:20.16f}"+"  "+f"{poscar1.lc[i][2]:20.16f}"+"\n")
#f2.write("volume "+f"{volume:24.16f}"+"\nend latticevecs\n\nbegin coordinates\n")
#k=0
#for i in range(poscar1.Ntype) :
#    f2.write("newtype "+poscar1.atomtype[i]+"\n")
#    for j in range(poscar1.Naint[i]) :
#        f2.write("coord  "+f"{poscar1.ap[k][0]:20.16f}"+"  "+f"{poscar1.ap[k][1]:20.16f}"+"  "+f"{poscar1.ap[k][2]:20.16f}"+"\n")
#        k=k+1
#f2.write("end coordinates\n\n")
#f2.close()

poscar1.filewrite_prt("prt.st")
