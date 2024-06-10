#!/usr/bin/env python3

# This script reads the afm simulation results and makes the plot

# The matrix in this script are mostly M[ny][nx] or M[nz][ny][nx] to match the size of image

import sys
import os
from commons import load_constant, load_palette
from classes import POSCAR
from v3math import v3pvm
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
import numpy as np

bohr=load_constant("bohr")
rydberg=load_constant("rydberg")
Ha=rydberg*2.0
electron=load_constant("electron")
angstrom=load_constant("angstrom")

latom=False
if "atom" in sys.argv :
    latom=True
    sys.argv.remove("atom")
lbohr=False
if "bohr" in sys.argv :
    lbohr=True
    sys.argv.remove("bohr")
ltilt=False
if "tilt" in sys.argv :
    ltilt=True
    sys.argv.remove("tilt")
icenter=1
for word in sys.argv :
    if word.isnumeric() :
        icenter=int(word)-1
        sys.argv.remove(word)

# ==================== read the input file ====================
x_spacing = 0.6
y_spacing = 0.6
z_spacing = 0.3
z_range = [5.7,6.3]
parallel = 1
k_spring = 0.8 # k in N/m
k_spring = k_spring * angstrom**2/electron # convert spring constant to eV/A^2

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
    elif word[0]=="parallel" :
        parallel=int(word[1])
    elif word[0]=="k_spring" :
        k_spring=float(word[1])
    elif word[0] in ["fdet","boundary"] :
        pass
    else :
        print("Warning: keyword \""+word[0]+"\" is not defined.")


x_spacing = x_spacing * bohr
y_spacing = y_spacing * bohr
z_spacing = z_spacing * bohr
x_range = [ x_range[0] * bohr, x_range[1] * bohr ]
y_range = [ y_range[0] * bohr, y_range[1] * bohr ]
z_range = [ z_range[0] * bohr, z_range[1] * bohr ]

# ==================== prepare the x and y coordinates ====================

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
z=z_range[0]
nz=0
zlist=[]
while (z<z_range[1]+1e-6) :
    zlist.append(z)
    z+=z_spacing
    nz+=1

f3=open("steps.dat","r")
line=f3.readlines()
f3.close()
steplist=[]
for l in line :
    word=l.split()
    if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
        continue
    steplist.append(int(word[0]))

lxincrease=True
k=0
ip=0
# movelist here is the index (ix,iy) of each point
movelist=[]
for iy in range(ny) :
    if lxincrease :
        for ix in range(nx) :
            if k==0 :
                movelist.append([])
            movelist[-1].append([ix,iy])
            k+=1
            if k==steplist[ip] :
                ip+=1
                k=0
    else :
        for ix in range(nx-1,-1,-1) :
            if k==0 :
                movelist.append([])
            movelist[-1].append([ix,iy])
            k+=1
            if k==steplist[ip] :
                ip+=1
                k=0
    lxincrease=not lxincrease

# ==================== calculate or read toten ====================

files = os.listdir()
# if the toten.dat exist, read it, if not, read from calculations outputs
if "toten.dat" in files :
    f5=open("toten.dat","r")
    line=f5.readlines()
    f5.close()
    il=1
    toten=[]
    for iz in range(nz) :
        toten1=[]
        for iy in range(ny) :
            toten0=[]
            for ix in range(nx) :
                toten0.append(float(line[il].split()[3]))
                il+=1
            toten1.append(toten0)
        toten.append(toten1)

else:
    # initialize toten matrix
    # toten[nz][ny][nx] in Ry
    toten=[]
    for iz in range(nz) :
        toten0=[]
        for iy in range(ny) :
            toten0.append([0.0]*nx)
        toten.append(toten0)

    for iz in range(nz) :
        for ip in range(parallel) :
            istep=-1
            filename="seq_"+str(iz+1)+"_"+str(ip+1)+"/parsec.out"
            f1=open(filename,"r")
            line=f1.readlines()
            f1.close()
            
            for l in line: 
                word=l.split()
                if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
                    continue
                if len(word)>=2 and word[0]=="Starting" and word[1]=="SCF..." :
                    istep+=1
                if len(word)>=5 and word[0]=="Total" and word[1]=="Energy" and word[2]=="=" :
                    toten[iz][movelist[ip][istep][1]][movelist[ip][istep][0]]=float(word[3])*rydberg # convert Ry to eV
    
    # write the toten file
    f5=open("toten.dat","w")
    f5.write("#ix iy iz toten(eV)\n")
    for iz in range(nz) :
        for iy in range(ny) :
            for ix in range(nx) :
                f5.write(str(ix+1)+" "+str(iy+1)+" "+str(iz+1)+" "+str(toten[iz][iy][ix])+"\n")
    f5.close()

# ==================== caclulate forces for tilt correction ====================

if ltilt :
    fx=[] # fx[ny][nx]
    fy=[] # fy[ny][nx]
    for iy in range(ny) :
        fx0=[0.0]*nx
        fx.append(fx0)
        fy0=[0.0]*nx
        fy.append(fy0)
    for iy in range(ny) :
        for ix in range(nx) :
            if ix!=0 and ix!=nx-1 :
                fx[iy][ix]=(toten[icenter][iy][ix-1]-toten[icenter][iy][ix+1])*0.5/x_spacing
            #elif ix==0 :
            #    fx[iy][ix]=(toten[1][iy][ix]-toten[1][iy][ix+1])/x_spacing
            #elif ix==nx-1 :
            #    fx[iy][ix]=(toten[1][iy][ix-1]-toten[1][iy][ix])/x_spacing
            if iy!=0 and iy!=ny-1 :
                fy[iy][ix]=(toten[icenter][iy-1][ix]-toten[icenter][iy+1][ix])*0.5/y_spacing
            #elif iy==0 :
            #    fy[iy][ix]=(toten[1][iy][ix]-toten[1][iy+1][ix])/y_spacing
            #elif iy==ny-1 :
            #    fy[iy][ix]=(toten[1][iy-1][ix]-toten[1][iy][ix])/y_spacing

    xy_new=[] # xy_new[ny][nx][2]
    for iy in range(ny) :
        xy_new0=[]
        for ix in range(nx) :
            xy_new0.append([ylist[iy] + fy[iy][ix] / k_spring, xlist[ix] + fx[iy][ix] / k_spring])
        xy_new.append(xy_new0)

    # create 2d map from toten
    toten_2d=[]
    for iz in range(nz) :
        toten_2d0 = scipy.interpolate.RectBivariateSpline(ylist, xlist, toten[iz])
        toten_2d.append(toten_2d0)

# ==================== caclulate kts ====================

kts=[]
if ltilt :
    for iy in range(ny) :
        kts0=[]
        for ix in range(nx) :
            y_new=xy_new[iy][ix][0]
            x_new=xy_new[iy][ix][1]
            #kts1=(0.25*toten_2d[icenter-1](y_new,x_new)[0,0]-0.5*toten_2d[icenter](y_new,x_new)[0,0]+0.25*toten_2d[icenter+1](y_new,x_new)[0,0])/z_spacing**2
            kts1=(toten_2d[icenter-1](y_new,x_new)[0,0]-2*toten_2d[icenter](y_new,x_new)[0,0]+toten_2d[icenter+1](y_new,x_new)[0,0])/z_spacing**2
            if lbohr :
                # convert k_ts from eV/A^2 to Ha/a0^2
                kts1 = kts1 * bohr**2/Ha
            kts0.append(kts1)
        kts.append(kts0)
else :
    for iy in range(ny) :
        kts0=[]
        for ix in range(nx) :
            #kts1=(0.25*toten[icenter-1][iy][ix]-0.5*toten[icenter][iy][ix]+0.25*toten[icenter+1][iy][ix])/z_spacing**2
            kts1=(toten[icenter-1][iy][ix]-2*toten[icenter][iy][ix]+toten[icenter+1][iy][ix])/z_spacing**2
            if lbohr :
                # convert k_ts from eV/A^2 to Ha/a0^2
                kts1 = kts1 * bohr**2/Ha
            kts0.append(kts1)
        kts.append(kts0)

# ==================== construct the atomic structure ====================
palette=load_palette()
if lbohr :
    funit=bohr
else : 
    funit=1.0
if latom :
    poscar1=POSCAR()
    poscar1.fileread_parsec(filename="sample.parsec_st.dat")
    poscar1.load_atomcolor()
    atom=[]
    color=[]
    ap=[]
    ia=0
    for itype in range(poscar1.Ntype) :
        for iatom in range(poscar1.Naint[itype]) :
            ap0=poscar1.ap[ia]
            for ix1 in [-1,0] :
                for ix2 in [-1,0] :
                    ap.append(v3pvm([ap0[0]+ix1,ap0[1]+ix2,ap0[2]],poscar1.lc))
                    atom.append(poscar1.atomtype[itype])
            ia+=1
    
    zmax=-1e6
    for ia in range(len(ap)) :
        zmax=max(zmax,ap[ia][2])

    atomx=[]
    atomy=[]
    atomcolor=[]
    edgecolor=[]
    for ia in range(len(ap)) :
        if ap[ia][2] > zmax-1.0 :
            atomx.append(ap[ia][0]/funit)
            atomy.append(ap[ia][1]/funit)
            atomcolor.append(palette[poscar1.atomcolor[atom[ia]]])
            if poscar1.atomcolor[atom[ia]].startswith("dark") or poscar1.atomcolor[atom[ia]]=="black" :
                edgecolor.append(palette["white"])
            else :
                edgecolor.append(palette["black"])

# ==================== creating the plot ====================
mpl.rcParams["font.sans-serif"].insert(0,"Noto Sans")
mpl.rcParams.update({'font.size': 14})
mpl.rcParams.update({'mathtext.default': 'regular'})

fig=plt.figure(figsize=(5,3.75))
gs0=fig.add_gridspec(1,2,wspace=0.02,hspace=0.00,left=0.14,right=0.80,top=0.95,bottom=0.15,width_ratios=[0.6,0.04])
[ax0,ax1]=gs0.subplots()

im_extent=[x_range[0]-x_spacing*0.5,x_range[1]+x_spacing*0.5,y_range[0]-y_spacing*0.5,y_range[1]+y_spacing*0.5]
for ic in range(len(im_extent)) :
    im_extent[ic]=im_extent[ic]/funit
im = ax0.imshow(kts, interpolation='bicubic', cmap="YlOrBr_r", origin="lower", extent=im_extent, aspect='equal', zorder=1)

if latom :
    ax0.scatter(atomx,atomy,c=atomcolor,s=12,edgecolors=edgecolor,linewidths=1,zorder=3)
ax0.set_xlim([x_range[0]/funit,x_range[1]/funit])
ax0.set_ylim([y_range[0]/funit,y_range[1]/funit])
if lbohr :
    ax0.set_xlabel("$\mathit{x}\ (Bohr)$",color=palette["black"])
    ax0.set_ylabel("$\mathit{y}\ (Bohr)$",color=palette["black"])
else :
    ax0.set_xlabel("$\mathit{x}\ (Å)$",color=palette["black"])
    ax0.set_ylabel("$\mathit{y}\ (Å)$",color=palette["black"])

cb=fig.colorbar(im, cax=ax1, orientation='vertical')
cb.outline.set_linewidth(1)
cb.outline.set_color(palette["black"])
if lbohr :
    ax1.set_ylabel("$\mathit{k}_{ts}\ (a.u.)$",color=palette["black"])
else :
    ax1.set_ylabel("$\mathit{k}_{ts}\ (eV/Å^2)$",color=palette["black"])

ax0.tick_params(axis="x", bottom=True, right=False, direction="in", color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
ax0.tick_params(axis="y", left=True, right=False, direction="in", color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
ax1.tick_params(axis="x", bottom=False, right=False, direction="in", color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
ax1.tick_params(axis="y", left=False, right=True, direction="in", color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
for edge in ["bottom", "top", "left", "right"] :
    ax0.spines[edge].set_color(palette["black"])
    ax0.spines[edge].set_linewidth(1)
    ax0.spines[edge].set_zorder(4)

filename="afm"
if ltilt :
    filename+="_tilt"
if latom :
    filename+="_atom"
if lbohr :
    filename+="_bohr"
filename+="_"+str(icenter+1)
filename+=".png"
fig.savefig(filename,dpi=1200)

