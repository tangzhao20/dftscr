#!/usr/bin/env python3

# This script works with afm.sh together to prepare the inputs for AFM simulation

from classes import POSCAR
from commons import load_constant, load_atom_dic
import sys
import os
import math

bohr=load_constant("bohr")
lvasp=False
if "vasp" in sys.argv :
    lvasp=True
    sys.argv.remove("vasp")

#================================================================
# read the input file
x_spacing=0.6
y_spacing=0.6
z_spacing=0.3
z_range=[5.7,6.3]
boundary=-1e0
lfdet=False

f1=open("afm.in","r")
line=f1.readlines()
f1.close()

for l in line :
    word=l.split()
    if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
        continue
    if word[0]=="x_range" :
        x_range=[float(word[1]),float(word[2])]
    elif word[0]=="y_range" :
        y_range=[float(word[1]),float(word[2])]
    elif word[0]=="z_range" :
        z_range=[float(word[1]),float(word[2])]
    elif word[0]=="x_spacing" :
        x_spacing=float(word[1])
    elif word[0]=="y_spacing" :
        y_spacing=float(word[1])
    elif word[0]=="z_spacing" :
        z_spacing=float(word[1])
    elif word[0]=="boundary" :
        boundary=float(word[1])
        if len(word)>2 and word[2][0] in ["a","A"] :
            boundary=boundary/bohr
    elif word[0]=="parallel" :
        parallel=int(word[1])
    elif word[0]=="fdet" :
        lfdet=True
        if len(word)>1 and word[1].lower() in ["false",".false."] :
            lfdet=False
    else :
        print("Warning: keyword "+word[0]+" is not defined.")

x=x_range[0]
nx=0
xlist=[]
while (x<x_range[1]+1e-6) :
    xlist.append(x)
    x+=x_spacing
    nx+=1
x_range[1]=xlist[-1]
y=y_range[0]
ny=0
ylist=[]
while (y<y_range[1]+1e-6) :
    ylist.append(y)
    y+=y_spacing
    ny+=1
y_range[1]=ylist[-1]
z=z_range[0]
nz=0
zlist=[]
while (z<z_range[1]+1e-6) :
    zlist.append(z)
    z+=z_spacing
    nz+=1
z_range[1]=zlist[-1]

#================================================================
# Read the structure from tip.xyz and sample.parsec_st.dat
poscar1=POSCAR()
poscar1.fileread_xyz("tip.xyz")
poscar2=POSCAR()
poscar2.fileread_parsec("sample.parsec_st.dat")

if poscar2.Ndim!=0 and poscar2.Ndim!=2 :
    print("Error: Only Ndim = 0 or 2 are supported.")
    sys.exit()

#================================================================
# Read the core charge from pp files *_POTRE.DAT
zion1=[]
for i in range(poscar1.Ntype) :
    pp_name=poscar1.atomtype[i]+"_POTRE.DAT"
    f5=open(pp_name,"r")
    line=f5.readlines()
    f5.close()
    zion1.append(int(float(line[3].split()[5])+0.5))
zion2=[]
for i in range(poscar2.Ntype) :
    pp_name=poscar2.atomtype[i]+"_POTRE.DAT"
    f5=open(pp_name,"r")
    line=f5.readlines()
    f5.close()
    zion2.append(int(float(line[3].split()[5])+0.5))

# Calculate the maximum of Nb for in advance.
atomdir=load_atom_dic()
z_full_list=[2,10,18,36,54,86,118]
Nb1_max=0
for i in range(poscar1.Ntype) :
    # Nb_atom=(z_full-z_atom+z_pp)/2
    z_atom=atomdir[poscar1.atomtype[i]]
    for j in range(7):
        if z_full_list[j] >= z_atom :
            z_full=z_full_list[j]
            break
    z_pp=zion1[i]
    Nb1_max+=(z_full-z_atom+z_pp)*poscar1.Naint[i]
Nb1_max=Nb1_max//2
Nb2_max=0
for i in range(poscar2.Ntype) :
    # Nb_atom=(z_full-z_atom+z_pp)/2
    z_atom=atomdir[poscar2.atomtype[i]]
    for j in range(7):
        if z_full_list[j] >= z_atom :
            z_full=z_full_list[j]
            break
    z_pp=zion2[i]
    Nb2_max+=(z_full-z_atom+z_pp)*poscar2.Naint[i]
Nb2_max=Nb2_max//2
# Calculate the prefered Nb of sample.
Nb2=0
for i in range(poscar2.Ntype) :
    Nb2+=zion2[i]*poscar2.Naint[i]
# Modify the following line if needed.
Nb2=int(Nb2/2*1.1+10)
Nb2=min(Nb2,Nb2_max)

# calculate the apc, this should be done before the loop
apc1=poscar1.cartesian(shift=[-0.5,-0.5,-0.5],factor=1.0/bohr)
if poscar2.Ndim==0 :
    shift=[-0.5,-0.5,-0.5]
elif poscar2.Ndim==2 :
    shift=[0.0,0.0,-0.5]
apc2=poscar2.cartesian(shift=shift,factor=1.0/bohr)

# move the tip atom to (0,0,0)
zmin=1e6
for ia in range(poscar1.Natom) :
    if zmin>apc1[ia][2] :
        zmin=apc1[ia][2]
        itip=ia
apmin=apc1[itip].copy()
# write poscar1.ap in Cartesian coordinate in Bohr
for ia in range(poscar1.Natom) :
    for ix in range(3) :
        apc1[ia][ix]=apc1[ia][ix]-apmin[ix]

# move the sample top to z=0
zmax=-1e6
for ia in range(poscar2.Natom) :
    if apc2[ia][2]>zmax :
        zmax=apc2[ia][2]

# instead of the very top atom, we select the medium of atoms with in 1 A under the top
z_toplist=[]
for ia in range(poscar2.Natom) :
    if apc2[ia][2]>zmax-1.0/bohr :
        z_toplist.append(apc2[ia][2])
z_toplist.sort()
ia_mid=len(z_toplist) // 2
zmax= (z_toplist[ia_mid] + z_toplist[~ia_mid]) / 2

for ia in range(poscar2.Natom) :
    apc2[ia][2]=apc2[ia][2]-zmax

#================================================================
# Write the manual.*.dat, which is splitted by the `parallel` parameter.

# calculate the step numbers. Need to minus 1 in the parsec.in file 
steplist=[(nx*ny)//parallel]*parallel
for ip in range((nx*ny)%parallel) :
    steplist[ip]+=1
f3=open("steps.dat","w")
for ip in range(parallel) :
    f3.write(str(steplist[ip])+"\n")
f3.close()
movelist=[]

lxincrease=True
k=0
ip=0
for iy in range(ny) :
    if lxincrease :
        for ix in range(nx) :
            if k==0 :
                movelist.append([])
            movelist[-1].append([xlist[ix],ylist[iy]])
            k+=1
            if k==steplist[ip] :
                ip+=1
                k=0
    else :
        for ix in range(nx-1,-1,-1) :
            if k==0 :
                movelist.append([])
            movelist[-1].append([xlist[ix],ylist[iy]])
            k+=1
            if k==steplist[ip] :
                ip+=1
                k=0
    lxincrease=not lxincrease

#================================================================
# If the boundary is not set, calculate it
if boundary<0 :
    if poscar2.Ndim==0 :
        for a in apc1 :
            for ix in range(2) :
                for iy in range(2) :
                    boundary=max(boundary,(a[0]+x_range[ix])**2+(a[1]+y_range[iy])**2+(a[2]+z_range[1])**2)
        for a in apc2 :
            boundary=max(boundary,a[0]**2+a[1]**2+a[2]**2)
        boundary=boundary**0.5+10.0

    elif poscar2.Ndim==2 :
        boundary_min=1e6
        boundary_max=-1e6
        for a in apc1 :
            boundary_max=max(boundary_max,a[2]+z_range[1])
        for a in apc2 :
            boundary_min=min(boundary_min,a[2])
        # Move the atoms to centralize the sample-tip structure
        z_move=0.5*(boundary_max+boundary_min)
        for ia in range(poscar1.Natom) :
            apc1[ia][2]-=z_move
        for ia in range(poscar2.Natom) :
            apc2[ia][2]-=z_move
        boundary=0.5*(boundary_max-boundary_min)+10.0
else : # Calculate the z_move, make the bottom at 10 bohr from boundary
    if poscar2.Ndim==2 :
        boundary_min=1e6
        for a in apc2 :
            boundary_min=min(boundary_min,a[2])
        z_move=boundary_min-(-boundary+10.0)
        for ia in range(poscar1.Natom) :
            apc1[ia][2]-=z_move
        for ia in range(poscar2.Natom) :
            apc2[ia][2]-=z_move

# if FDET, we do a calculation for a sample potential file
if lfdet :
    filename_parsec="parsec_st_spot.dat"
    f2=open(filename_parsec,"w")
    f2.write("#---------output from afm.py----------\n")
    f2.write("States_Num "+str(Nb2)+"\n\n")
    if poscar2.Ndim==0:
        pass
    if poscar2.Ndim==2 :
        f2.write("Boundary_Conditions slab\n")
        f2.write("begin Cell_Shape\n")
        for ix1 in range(2) :
            for ix2 in range(3) :
                f2.write(f"{poscar2.lc[ix1][ix2]/bohr:18.12f}")
            f2.write("\n")
        f2.write("end Cell_Shape\n\n")

        f2.write("Kpoint_Method mp\n\n")
        f2.write("begin Monkhorst_Pack_Grid\n")
        for ix in range(2) :
            kgrid=math.floor(30.0/(poscar2.lc[ix][0]**2+poscar2.lc[ix][1]**2+poscar2.lc[ix][2]**2)**0.5)+1
            f2.write(f"  {kgrid:d}")
        f2.write("\nend Monkhorst_Pack_Grid\n\n")
        f2.write("begin Monkhorst_Pack_Shift\n")
        f2.write("0.0  0.0  0.0\n")
        f2.write("end Monkhorst_Pack_Shift\n\n")

    f2.write(f"Boundary_Sphere_Radius {boundary:.12g}\n\n")
    f2.write("Atom_Types_Num "+str(poscar2.Ntype)+"\n")
    f2.write("Coordinate_Unit Cartesian_Bohr\n\n")

    f2.write("#------------ begin sample -----------\n")
    k=0
    for i in range(poscar2.Ntype) :
        f2.write("Atom_Type: "+poscar2.atomtype[i]+"\n")
        f2.write("Local_Component: s\n")
        f2.write("begin Atom_Coord\n")
        for ia in range(poscar2.Naint[i]) :
            f2.write(f"{apc2[k][0]:18.12f}{apc2[k][1]:18.12f}{apc2[k][2]:18.12f}\n")
            k=k+1
        f2.write("end Atom_Coord\n\n")
    f2.write("#------------- end sample ------------\n\n")

    f2.close()

# main loop
for iz in range(nz) :
    for ip in range(parallel) :
        #================================================================
        # Write the structure to parsec_st_iz_ip.dat file
        filename_parsec="parsec_st_"+str(iz+1)+"_"+str(ip+1)+".dat"
        f2=open(filename_parsec,"w")
        f2.write("#---------output from afm.py----------\n")
        if lfdet :
            f2.write("States_Num "+str(Nb1_max)+"\n\n")
        else :
            f2.write("States_Num "+str(Nb1_max+Nb2)+"\n\n")
        if poscar2.Ndim==0:
            pass
        if poscar2.Ndim==2 :
            f2.write("Boundary_Conditions slab\n")
            f2.write("begin Cell_Shape\n")
            for ix1 in range(2) :
                for ix2 in range(3) :
                    f2.write(f"{poscar2.lc[ix1][ix2]/bohr:18.12f}")
                f2.write("\n")
            f2.write("end Cell_Shape\n\n")

            f2.write("Kpoint_Method mp\n\n")
            f2.write("begin Monkhorst_Pack_Grid\n")
            for ix in range(2) :
                kgrid=math.floor(30.0/(poscar2.lc[ix][0]**2+poscar2.lc[ix][1]**2+poscar2.lc[ix][2]**2)**0.5)+1
                f2.write(f"  {kgrid:d}")
            f2.write("\nend Monkhorst_Pack_Grid\n\n")
            f2.write("begin Monkhorst_Pack_Shift\n")
            f2.write("0.0  0.0  0.0\n")
            f2.write("end Monkhorst_Pack_Shift\n\n")

        f2.write(f"Boundary_Sphere_Radius {boundary:.12g}\n\n")
        if lfdet :
            f2.write("Atom_Types_Num "+str(poscar1.Ntype)+"\n")
        else :
            f2.write("Atom_Types_Num "+str(poscar1.Ntype+poscar2.Ntype)+"\n")
        f2.write("Coordinate_Unit Cartesian_Bohr\n\n")
        
        f2.write("#------------- begin tip -------------\n")
        k=0
        for i in range(poscar1.Ntype) :
            f2.write("Atom_Type: "+poscar1.atomtype[i]+"\n")
            f2.write("Local_Component: s\n")  # Assume the local_component is s here. A dictionary of {element: local_component} should be added
            f2.write("begin Atom_Coord\n")
            for ia in range(poscar1.Naint[i]) :
                f2.write(f"{apc1[k][0]+movelist[ip][0][0]:18.12f}{apc1[k][1]+movelist[ip][0][1]:18.12f}{apc1[k][2]+zlist[iz]:18.12f}\n")
                k=k+1
            f2.write("end Atom_Coord\n\n")
        f2.write("#-------------- end tip --------------\n\n")
        
        f2.write("#------------ begin sample -----------\n")
        if lfdet :
            f2.write("Add_Point_Charges: .TRUE.\n")
            f2.write("Point_Typ_Num: "+str(poscar2.Ntype)+"\n\n")

        k=0
        for i in range(poscar2.Ntype) :
            if lfdet:
                f2.write("Pt_Chg: "+str(zion2[i])+"\n")
                f2.write("begin Point_Coord\n")
            else:
                f2.write("Atom_Type: "+poscar2.atomtype[i]+"\n")
                f2.write("Local_Component: s\n")
                f2.write("begin Atom_Coord\n")

            for ia in range(poscar2.Naint[i]) :
                f2.write(f"{apc2[k][0]:18.12f}{apc2[k][1]:18.12f}{apc2[k][2]:18.12f}\n")
                k=k+1

            if lfdet:
                f2.write("end Point_Coord\n\n")
            else :
                f2.write("end Atom_Coord\n\n")
        f2.write("#------------- end sample ------------\n\n")
        
        f2.close()

        # convert the structure to vasp format
        if lvasp :
            os.system("posconvert.py parsec vasp "+filename_parsec)
            os.system("mv POSCAR.new "+str(iz+1)+"_"+str(ip+1)+".vasp")

        # write the manual_*.dat file
        f4=open("manual_"+str(iz+1)+"_"+str(ip+1)+".dat","w")
        for istep in range(1,steplist[ip]) :
            for ia in range(poscar1.Natom) :
                f4.write(f"{apc1[ia][0]+movelist[ip][istep][0]:18.12f}{apc1[ia][1]+movelist[ip][istep][1]:18.12f}{apc1[ia][2]+zlist[iz]:18.12f}\n")
            if not lfdet :
                for ia in range(poscar2.Natom) :
                    f4.write(f"{apc2[ia][0]:18.12f}{apc2[ia][1]:18.12f}{apc2[ia][2]:18.12f}\n")
            f4.write("\n")
        f4.close()
 
