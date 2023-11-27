#!/usr/bin/env python3

# from the conventional cell of silicon
# 1. create a supercell, whose size is determined by raidus
# 2. select the atoms inside of the sphere
# 3. add hydrogen atoms

import sys
from classes import POSCAR

radius=float(sys.argv[1])
print(" raidus = "+str(radius)+" ")

si_a0=[[0.0,0.0,0.0],[0.0,0.5,0.5],[0.5,0.0,0.5],[0.5,0.5,0.0]] # the site-A silicon atoms at (0,0,0)
si_b0=[[0.25,0.25,0.25],[0.25,0.75,0.75],[0.75,0.25,0.75],[0.75,0.75,0.25]] # the site-B silicon atoms at (1/4,1/4,1/4)
#si_a0=[[-0.125,-0.125,-0.125],[-0.125,0.375,0.375],[0.375,-0.125,0.375],[0.375,0.375,-0.125]] # the site-A silicon atoms at (-1/8,-1/8,-1/8)
#si_b0=[[0.125,0.125,0.125],[0.125,0.625,0.625],[0.625,0.125,0.625],[0.625,0.625,0.125]] # the site-B silicon atoms at (1/8,1/8,1/8)
H=[] # hydrogen

#lc=5.431020511 # Ang, lenc
lc=5.43
lc_4=lc/4
si_h_4=1.46/3**0.5

Ncell=int(radius/lc)+1

# Set up the Si supercell
si_a2=[]
si_b2=[]
for i0 in range(-Ncell,Ncell+1) :
    for i1 in range(-Ncell,Ncell+1) :
        for i2 in range(-Ncell,Ncell+1) :
            for j in range(4) :
                si_a2.append([(si_a0[j][0]+i0)*lc,(si_a0[j][1]+i1)*lc,(si_a0[j][2]+i2)*lc])
                si_b2.append([(si_b0[j][0]+i0)*lc,(si_b0[j][1]+i1)*lc,(si_b0[j][2]+i2)*lc])

# Remove Si atoms from outside the sphere
si_a=[]
for si in si_a2 :
    if si[0]**2+si[1]**2+si[2]**2 <= radius**2 :
        si_a.append(si)
si_b=[]
for si in si_b2 :
    if si[0]**2+si[1]**2+si[2]**2 <= radius**2 :
        si_b.append(si)

# setup Si-Si bonds
# -1: no atom, i: the ith atom in the other list
si_a_bond=[] # (+++,+--,-+-,--+)
for ia in range(len(si_a)) :
    si_a_bond.append([-1]*4)
si_b_bond=[] # (---,-++,+-+,++-)
for ia in range(len(si_b)) :
    si_b_bond.append([-1]*4)

for i1 in range(len(si_a)) :
    for i2 in range(len(si_b)) :
        if max(abs(si_b[i2][0]-si_a[i1][0]-lc_4), abs(si_b[i2][1]-si_a[i1][1]-lc_4), abs(si_b[i2][2]-si_a[i1][2]-lc_4))<=1e-6 :
            si_a_bond[i1][0]=i2
            si_b_bond[i2][0]=i1
        elif max(abs(si_b[i2][0]-si_a[i1][0]-lc_4), abs(si_b[i2][1]-si_a[i1][1]+lc_4), abs(si_b[i2][2]-si_a[i1][2]+lc_4))<=1e-6 :
            si_a_bond[i1][1]=i2
            si_b_bond[i2][1]=i1
        elif max(abs(si_b[i2][0]-si_a[i1][0]+lc_4), abs(si_b[i2][1]-si_a[i1][1]-lc_4), abs(si_b[i2][2]-si_a[i1][2]+lc_4))<=1e-6 :
            si_a_bond[i1][2]=i2
            si_b_bond[i2][2]=i1
        elif max(abs(si_b[i2][0]-si_a[i1][0]+lc_4), abs(si_b[i2][1]-si_a[i1][1]+lc_4), abs(si_b[i2][2]-si_a[i1][2]-lc_4))<=1e-6 :
            si_a_bond[i1][3]=i2
            si_b_bond[i2][3]=i1

# Remove Si atoms with only 1 Si-Si bond
while True :
    fsingle=False
    nneib=[0]*len(si_a)
    inew=[-1]*len(si_a)
    ilog=0
    for i1 in range(len(si_a)) :
        for j in range(4) :
            if si_a_bond[i1][j]!=-1 :
                nneib[i1]+=1
        if nneib[i1]<=1 :
            fsingle=True
        else :
            inew[i1]=ilog
            ilog+=1
    for i1 in range(len(si_b)) :
        for j in range(4) :
            if si_b_bond[i1][j]>=0 :
                si_b_bond[i1][j]=inew[si_b_bond[i1][j]]
    si_a=[si_a[i1] for i1 in range(len(si_a)) if nneib[i1]>1 ]
    si_a_bond=[si_a_bond[i1] for i1 in range(len(si_a_bond)) if nneib[i1]>1 ]

    nneib=[0]*len(si_b)
    inew=[-1]*len(si_b)
    ilog=0
    for i1 in range(len(si_b)) :
        for j in range(4) :
            if si_b_bond[i1][j]!=-1 :
                nneib[i1]+=1
        if nneib[i1]<=1 :
            fsingle=True
        else :
            inew[i1]=ilog
            ilog+=1
    for i1 in range(len(si_a)) :
        for j in range(4) :
            if si_a_bond[i1][j]>=0 :
                si_a_bond[i1][j]=inew[si_a_bond[i1][j]]
    si_b=[si_b[i1] for i1 in range(len(si_b)) if nneib[i1]>1 ]
    si_b_bond=[si_b_bond[i1] for i1 in range(len(si_b_bond)) if nneib[i1]>1 ]
    if fsingle==False :
        break

# Add H atoms
for ia in range(len(si_a)) :
    if si_a_bond[ia][0]==-1 :
        H.append([si_a[ia][0]+si_h_4, si_a[ia][1]+si_h_4, si_a[ia][2]+si_h_4])
    if si_a_bond[ia][1]==-1 :
        H.append([si_a[ia][0]+si_h_4, si_a[ia][1]-si_h_4, si_a[ia][2]-si_h_4])
    if si_a_bond[ia][2]==-1 :
        H.append([si_a[ia][0]-si_h_4, si_a[ia][1]+si_h_4, si_a[ia][2]-si_h_4])
    if si_a_bond[ia][3]==-1 :
        H.append([si_a[ia][0]-si_h_4, si_a[ia][1]-si_h_4, si_a[ia][2]+si_h_4])
for ia in range(len(si_b)) :
    if si_b_bond[ia][0]==-1 :
        H.append([si_b[ia][0]-si_h_4, si_b[ia][1]-si_h_4, si_b[ia][2]-si_h_4])
    if si_b_bond[ia][1]==-1 :
        H.append([si_b[ia][0]-si_h_4, si_b[ia][1]+si_h_4, si_b[ia][2]+si_h_4])
    if si_b_bond[ia][2]==-1 :
        H.append([si_b[ia][0]+si_h_4, si_b[ia][1]-si_h_4, si_b[ia][2]+si_h_4])
    if si_b_bond[ia][3]==-1 :
        H.append([si_b[ia][0]+si_h_4, si_b[ia][1]+si_h_4, si_b[ia][2]-si_h_4])

# write to file
poscar1=POSCAR(empty=True)
lc_big=2*radius+10
poscar1.lc=[[lc_big,0,0],[0,lc_big,0],[0,0,lc_big]]
poscar1.Natom=len(si_a)+len(si_b)+len(H)
poscar1.Ntype=2
poscar1.atomtype=["Si","H"]
poscar1.Naint=[len(si_a)+len(si_b),len(H)]
for a in si_a :
    poscar1.ap.append([a[0]/lc_big+0.5,a[1]/lc_big+0.5,a[2]/lc_big+0.5])
for a in si_b :
    poscar1.ap.append([a[0]/lc_big+0.5,a[1]/lc_big+0.5,a[2]/lc_big+0.5])
for a in H :
    poscar1.ap.append([a[0]/lc_big+0.5,a[1]/lc_big+0.5,a[2]/lc_big+0.5])
poscar1.movetobox()
poscar1.filewrite("Si"+str(poscar1.Naint[0])+"H"+str(poscar1.Naint[1])+".vasp")
