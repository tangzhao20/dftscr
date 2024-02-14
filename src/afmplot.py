#!/usr/bin/env python3

# This script reads the afm simulation results and makes the plot

import sys
import os
from commons import load_palette
from classes import POSCAR
from v3math import v3pvm
import matplotlib as mpl
import matplotlib.pyplot as plt
#import matplotlib.cm as cm

latom=False
if "atom" in sys.argv :
    latom=True
    sys.argv.remove("atom")

# ==================== read the input file ====================
x_spacing=0.6
y_spacing=0.6
z_sampling=[5.7,6.0,6.3]
z_spacing=0.3
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

files = os.listdir()
fkts=False
# if the kts.dat exist, read it, if not, read from calculations outputs
if "kts.dat" in files :
    f4=open("kts.dat","r")
    line=f4.readlines()
    f4.close()
    il=0
    kts=[]
    for ix in range(nx) :
        kts0=[]
        for iy in range(ny) :
            kts0.append(float(line[il].split()[2]))
            il+=1
        kts.append(kts0)

else:
    # initialize toten matrix
    toten=[]
    for iz in range(3) :
        toten0=[]
        for iy in range(ny) :
            toten0.append([0.0]*nx)
        toten.append(toten0)

    for iz in range(3) :
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
                    toten[iz][movelist[ip][istep][1]][movelist[ip][istep][0]]=float(word[3]) # in Ry

    kts=[]
    for ix in range(nx) :
        kts0=[]
        for iy in range(ny) :
            # why this is not (2E1-E0-E2)/z^2?
            kts1=(0.25*toten[0][ix][iy]-0.5*toten[1][ix][iy]+0.25*toten[2][ix][iy])/z_spacing**2*1000
            kts0.append(kts1)
        kts.append(kts0)
    kts.reverse() # Need to check: imshow needs the matrix with reversed rows?

    f4=open("kts.dat","w")
    for ix in range(nx) :
        for iy in range(ny) :
            f4.write(str(ix+1)+" "+str(iy+1)+" "+str(kts[ix][iy])+"\n")
    f4.close()

palette=load_palette() 
# ==================== construct the atomic structure ====================
if latom :
    poscar1=POSCAR(empty=True)
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


    bohr=0.529177249
    atomx=[]
    atomy=[]
    atomcolor=[]
    edgecolor=[]
    for ia in range(len(ap)) :
        if ap[ia][2] > zmax-1.0 :
            atomx.append(ap[ia][0]/bohr)
            atomy.append(ap[ia][1]/bohr)
            atomcolor.append(palette[poscar1.atomcolor[atom[ia]]])
            if poscar1.atomcolor[atom[ia]].startswith("dark") or poscar1.atomcolor[atom[ia]]=="black" :
                edgecolor.append(palette["white"])
            else :
                edgecolor.append(palette["black"])

mpl.rcParams["font.sans-serif"].insert(0,"Noto Sans")
mpl.rcParams.update({'font.size': 14})
mpl.rcParams.update({'mathtext.default': 'regular'})

fig=plt.figure(figsize=(5,3.75))
gs0=fig.add_gridspec(1,2,wspace=0.02,hspace=0.00,left=0.14,right=0.80,top=0.95,bottom=0.15,width_ratios=[0.6,0.04])
[ax0,ax1]=gs0.subplots()

#im=ax0.imshow(kts,interpolation='gaussian',cmap=cm.gray,extent=[x_range[0],x_range[1],y_range[0],y_range[1]],aspect='equal',zorder=1)
im=ax0.imshow(kts,interpolation='bicubic',cmap="YlOrBr_r",extent=[x_range[0],x_range[1],y_range[0],y_range[1]],aspect='equal',zorder=1)

if latom :
    ax0.scatter(atomx,atomy,c=atomcolor,s=12,edgecolors=edgecolor,linewidths=1,zorder=3)
ax0.set_xlim([x_range[0],x_range[1]])
ax0.set_ylim([y_range[0],y_range[1]])
ax0.set_xlabel("$\mathit{x}\ (Bohr)$",color=palette["black"])
ax0.set_ylabel("$\mathit{y}\ (Bohr)$",color=palette["black"])

cb=fig.colorbar(im, cax=ax1, orientation='vertical')
cb.outline.set_linewidth(1)
cb.outline.set_color(palette["black"])
ax1.set_ylabel("$\mathit{k}_{ts}\ (10^{-3}\ a.u.)$",color=palette["black"])

ax0.tick_params(axis="x", bottom=True, right=False, direction="in", color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
ax0.tick_params(axis="y", left=True, right=False, direction="in", color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
ax1.tick_params(axis="x", bottom=False, right=False, direction="in", color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
ax1.tick_params(axis="y", left=False, right=True, direction="in", color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
for edge in ["bottom", "top", "left", "right"] :
    ax0.spines[edge].set_color(palette["black"])
    ax0.spines[edge].set_linewidth(1)
    ax0.spines[edge].set_zorder(4)

if latom :
    fig.savefig("afmatom.png",dpi=1200)
else :
    fig.savefig("afm.png",dpi=1200)


