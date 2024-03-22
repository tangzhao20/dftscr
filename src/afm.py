#!/usr/bin/env python3

# This script works with afm.sh together to prepare the inputs for AFM simulation

from classes import POSCAR
import sys
import os
import math

#================================================================
# read the input file
x_spacing=0.6
y_spacing=0.6
z_spacing=0.3
z_range=[5.7,6.3]
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
        z_range=[float(word[0]),float(word[1])]
    elif word[0]=="x_spacing" :
        x_spacing=float(word[1])
    elif word[0]=="y_spacing" :
        y_spacing=float(word[1])
    elif word[0]=="z_spacing" :
        z_spacing=float(word[1])
    elif word[0]=="parallel" :
        parallel=int(word[1])
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
# Read the structure from .xyz files
bohr=0.529177210903
poscar1=POSCAR(empty=True)
poscar1.fileread_xyz("tip.xyz")
poscar2=POSCAR(empty=True)
poscar2.fileread_parsec("sample.parsec_st.dat")

# calculate the apc, this should be done before the loop
if poscar2.Ndim==2 :
    shift=[0.0,0.0,0.5]
elif poscar2.Ndim==0 :
    shift=[0.5,0.5,0.5]
k=0
apc1=[]
for i in range(poscar1.Ntype) :
    for ia in range(poscar1.Naint[i]) :
        apc0=[0.0]*3
        for ix1 in range(3) :
            for ix2 in range(3) :
                apc0[ix2]+=(poscar1.ap[k][ix1]-shift[ix1])*poscar1.lc[ix1][ix2]
        for ix in range(3) :
            apc0[ix]=apc0[ix]/bohr
        apc1.append(apc0)
        k=k+1
k=0
apc2=[]
for i in range(poscar2.Ntype) :
    for ia in range(poscar2.Naint[i]) :
        apc0=[0.0]*3
        for ix1 in range(3) :
            for ix2 in range(3) :
                apc0[ix2]+=(poscar2.ap[k][ix1]-shift[ix1])*poscar2.lc[ix1][ix2]
        for ix in range(3) :
            apc0[ix]=apc0[ix]/bohr
        apc2.append(apc0)
        k=k+1

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
# Calculate the radius

if poscar2.Ndim==0:
    radius=0.0
    for a in apc1 :
        for ix in range(2) :
            for iy in range(2) :
                radius=max(radius,(a[0]+x_range[ix])**2+(a[1]+y_range[iy])**2+(a[2]+z_range[1])**2)
    for a in apc2 :
        radius=max(radius,a[0]**2+a[1]**2+a[2]**2)
    radius=radius**0.5+10.0

elif poscar2.Ndim==2 :
    radius_min=1e6
    radius_max=-1e6
    for a in apc1 :
        radius_max=max(radius_max,a[2]+z_range[1])
    for a in apc2 :
        radius_min=min(radius_min,a[2])
    # Move the atoms to centralize the sample-tip structure
    z_move=0.5*(radius_max+radius_min)
    for ia in range(poscar1.Natom) :
        apc1[ia][2]-=z_move
    for ia in range(poscar2.Natom) :
        apc2[ia][2]-=z_move
    radius=0.5*(radius_max-radius_min)+10.0

for iz in range(nz) :
    for ip in range(parallel) :
        #================================================================
        # Write the structure to parsec_st_iz_ip.dat file
        filename_parsec="parsec_st_"+str(iz+1)+"_"+str(ip+1)+".dat"
        f2=open(filename_parsec,"w")
        f2.write("#---------output from afm.py----------\n")
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
        else :
            print("Error: Only Ndim = 0 or 2 are supported.")
            sys.exit()
        
        f2.write(f"Boundary_Sphere_Radius {radius:.12g}\n\n")
        f2.write("Atom_Types_Num "+str(len(poscar1.atomtype)+len(poscar2.atomtype))+"\n")
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

        # convert the structure to vasp format
        if len(sys.argv)>1 and sys.argv[1]=="vasp" :
            os.system("posconvert.py parsec vasp "+filename_parsec)
            os.system("mv POSCAR.new "+str(iz+1)+"_"+str(ip+1)+".vasp")

        # write the manual_*.dat file
        f4=open("manual_"+str(iz+1)+"_"+str(ip+1)+".dat","w")
        for istep in range(1,steplist[ip]) :
            for ia in range(poscar1.Natom) :
                f4.write(f"{apc1[ia][0]+movelist[ip][istep][0]:18.12f}{apc1[ia][1]+movelist[ip][istep][1]:18.12f}{apc1[ia][2]+zlist[iz]:18.12f}\n")
            for ia in range(poscar2.Natom) :
                f4.write(f"{apc2[ia][0]:18.12f}{apc2[ia][1]:18.12f}{apc2[ia][2]:18.12f}\n")
            f4.write("\n")
        f4.close()
 
