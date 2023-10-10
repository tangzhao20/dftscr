#!/usr/bin/env python3

# This script works with afm.sh together to prepare the inputs for AFM simulation

from classes import POSCAR
import sys
import math

#================================================================
# read the input file
x_spacing=0.6
y_spacing=0.6
z_sampling=[5.7,6.0,6.3]
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
    elif word[0]=="x_spacing" :
        x_spacing=float(word[1])
    elif word[0]=="y_spacing" :
        y_spacing=float(word[1])
    elif word[0]=="z_sampling" :
        z_sampling=word[1:]
        for iz in range(range(z_sampling)) :
            try :
                z_sampling[iz]=float(z_sampling[iz])
            except :
                del zsampling[iz:]
                break
    elif word[0]=="parallel" :
        parallel=int(word[1])
    else :
        print("Warning: keyword "+word[0]+" is not defined.")

#================================================================
# Read the structure from .xyz files
bohr=0.529177210903
poscar1=POSCAR(empty=True)
poscar1.fileread_parsec("sample.parsec_st.dat")
poscar2=POSCAR(empty=True)
poscar2.fileread_xyz("tip.xyz")

# move the sample top to z=0
zmax=-1e6
for ia in range(poscar1.Natom) :
    zmax=max(poscar1.ap[ia][2],zmax)

# write poscar2.ap in Cartesian coordinate in Bohr
lc1=poscar1.lc[0][0]
for ia in range(poscar1.Natom) :
    poscar1.ap[ia][0]=(poscar1.ap[ia][0]-0.5)*lc1/bohr
    poscar1.ap[ia][1]=(poscar1.ap[ia][1]-0.5)*lc1/bohr
    poscar1.ap[ia][2]=(poscar1.ap[ia][2]-zmax)*lc1/bohr

# move the tip atom to (0,0,zz)
zmin=1e6
for ia in range(poscar2.Natom) :
    if zmin>poscar2.ap[ia][2] :
        zmin=poscar2.ap[ia][2]
        itip=ia
apmin=poscar2.ap[itip].copy()

# write poscar2.ap in Cartesian coordinate in Bohr
lc2=poscar2.lc[0][0]
for ia in range(poscar2.Natom) :
    poscar2.ap[ia][0]=(poscar2.ap[ia][0]-apmin[0])*lc2/bohr
    poscar2.ap[ia][1]=(poscar2.ap[ia][1]-apmin[1])*lc2/bohr
    poscar2.ap[ia][2]=(poscar2.ap[ia][2]-apmin[2])*lc2/bohr

#================================================================
# Write the structure to parsec_st.dat file
f2=open("parsec_st.dat","w")
f2.write("#---------output from afm.py----------\n")
radius=0.0
if poscar1.ndim==0:
    for a in poscar1.ap :
        radius=max(radius,a[0]**2+a[1]**2+a[2]**2)
    for a in poscar2.ap :
        radius=max(radius,(a[0]+x_range[0])**2+(a[1]+y_range[0])**2+(a[2]+max(z_sampling))**2)
    for a in poscar2.ap :
        radius=max(radius,(a[0]+x_range[1])**2+(a[1]+y_range[0])**2+(a[2]+max(z_sampling))**2)
    for a in poscar2.ap :
        radius=max(radius,(a[0]+x_range[0])**2+(a[1]+y_range[1])**2+(a[2]+max(z_sampling))**2)
    for a in poscar2.ap :
        radius=max(radius,(a[0]+x_range[1])**2+(a[1]+y_range[1])**2+(a[2]+max(z_sampling))**2)
    radius=radius**0.5+10.0
elif poscar1.ndim==2 :
    for a in poscar1.ap :
        radius=max(radius,abs(a[2]))
    for a in poscar2.ap :
        radius=max(radius,abs(a[2]+max(z_sampling)))
    radius=radius+10.0
    f2.write("Boundary_Conditions slab\n")
    f2.write("begin Cell_Shape\n")
    for ix1 in range(2) :
        for ix2 in range(3) :
            f2.write(f"{poscar1.lc[ix1][ix2]/bohr:18.12f}")
        f2.write("\n")
    f2.write("end Cell_Shape\n\n")

    f2.write("Kpoint_Method mp\n\n")
    f2.write("begin Monkhorst_Pack_Grid\n")
    for ix in range(2) :
        kgrid=math.floor(30.0/(poscar1.lc[ix][0]**2+poscar1.lc[ix][1]**2+poscar1.lc[ix][2]**2)**0.5)+1
        f2.write(f"  {kgrid:d}")
    f2.write("\nend Monkhorst_Pack_Grid\n\n")
    f2.write("begin Monkhorst_Pack_Shift\n")
    f2.write("0.0  0.0  0.0\n")
    f2.write("end Monkhorst_Pack_Shift\n\n")
else :
    print("Error: Only ndim = 0 or 2 are supported.")
    sys.exit()

f2.write(f"Boundary_Sphere_Radius {radius:.12g}\n\n")
f2.write("Atom_Types_Num "+str(len(poscar1.atomtype)+len(poscar2.atomtype))+"\n")
f2.write("Coordinate_Unit Cartesian_Bohr\n\n")

f2.write("#------------- begin tip -------------\n")
k=0
for i in range(poscar2.Ntype) :
    f2.write("Atom_Type: "+poscar2.atomtype[i]+"\n")
    f2.write("Local_Component: \n")
    f2.write("begin Atom_Coord\n")
    for ia in range(poscar2.Naint[i]) :
        f2.write(f"{poscar2.ap[k][0]:18.12f}{poscar2.ap[k][1]:18.12f}{poscar2.ap[k][2]+z_sampling[1]:18.12f}\n")
        k=k+1
    f2.write("end Atom_Coord\n\n")
f2.write("#-------------- end tip --------------\n\n")

f2.write("#------------ begin sample -----------\n")
k=0
for i in range(poscar1.Ntype) :
    f2.write("Atom_Type: "+poscar1.atomtype[i]+"\n")
    f2.write("Local_Component: \n")
    f2.write("begin Atom_Coord\n")
    for ia in range(poscar1.Naint[i]) :
        f2.write(f"{poscar1.ap[k][0]:18.12f}{poscar1.ap[k][1]:18.12f}{poscar1.ap[k][2]:18.12f}\n")
        k=k+1
    f2.write("end Atom_Coord\n\n")
f2.write("#------------- end sample ------------\n\n")

f2.close()

#================================================================
# Write the manual.*.dat, which is splitted by the `parallel` parameter.

x=x_range[0]
nx=0
xlist=[]
while (x<x_range[1]+1e-6) :
    xlist.append(x)
    x+=x_spacing
    nx+=1
y=y_range[0]
ny=0
ylist=[]
while (y<y_range[1]+1e-6) :
    ylist.append(y)
    y+=y_spacing
    ny+=1

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

for iz in range(3) :
    for ip in range(parallel) :
        f4=open("manual_"+str(iz+1)+"_"+str(ip+1)+".dat","w")
        for istep in range(steplist[ip]) :
            for a in poscar2.ap :
                f4.write(f"{a[0]+movelist[ip][istep][0]:18.12f}{a[1]+movelist[ip][istep][1]:18.12f}{a[2]+z_sampling[iz]:18.12f}\n")
            for a in poscar1.ap :
                f4.write(f"{a[0]:18.12f}{a[1]:18.12f}{a[2]:18.12f}\n")
            f4.write("\n")
        f4.close()
 
