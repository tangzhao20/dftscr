#!/usr/bin/env python3

# plot DOS

# python dos.py (v) package (E1) (E2)

# Input: DOSCAR

import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from classes import doscar
from commons import load_packagename, load_palette

fvertical=False
filename_out="dos.png"
if len(sys.argv)>1 and sys.argv[1] in ["v","vertical"] :
    fvertical=True
    filename_out="dos_v.png"
    del sys.argv[1]

if len(sys.argv)<=1 :
    print("python dos.py (v) package (Emin) (Emax)")
    sys.exit()
package=sys.argv[1]

packagename=load_packagename()

if package in packagename["vasp"]: 
    doscar0=doscar.DOSCAR()
elif package in packagename["qe"]:
    doscar0=doscar.DOSCAR(empty=True)
    # find a *.dos file
    files = os.listdir(".")
    for f in files:
        if f.endswith('.dos'):
            filename=f
            break
    doscar0.fileread_qe(filename)
else:
    print("Package \""+package+"\" is not supported yet.")
    print("python dos.py (v) package (Emin) (Emax)")
    sys.exit()

if len(sys.argv)>=4 :
    xmax=float(sys.argv[3])
    xmin=float(sys.argv[2])
elif len(sys.argv)==3 :
    xmax=float(sys.argv[2])
    xmin=-float(sys.argv[2])
else :
    xmax=5.0
    xmin=-5.0

doscar0.energyshift(doscar0.ef)

dosmax=0.0
for ie in range(doscar0.Nedos) :
    if doscar0.energy[ie]<xmin :
        continue
    if doscar0.energy[ie]>xmax :
        break
    if doscar0.Ns==1 :
        dosmax=max(dosmax,doscar0.dos[0][ie])
    else :
        dosmax=max(dosmax,doscar0.dos[0][ie],doscar0.dos[1][ie])

colpal=load_palette() #color palette: [blue, orange, gray, white, black]
mpl.rcParams["font.sans-serif"].insert(0,"Noto Sans")
mpl.rcParams.update({'font.size': 14})

if fvertical==False :
    fig=plt.figure(figsize=(5,3.75))
    if doscar0.Ns==1 :
        gs0=fig.add_gridspec(1,1,wspace=0.0,hspace=0.00,left=0.17,right=0.95,top=0.97, bottom=0.115)
        ax0=gs0.subplots()
    else :
        gs0=fig.add_gridspec(2,1,wspace=0.0,hspace=0.00,left=0.17,right=0.95,top=0.97, bottom=0.115)
        (ax0,ax1)=gs0.subplots()

    ax0.axvline(linewidth=1,color=colpal[2],zorder=0)
    ax0.plot(doscar0.energy,doscar0.dos[0],color=colpal[0],linewidth=1,zorder=3)
    ax0.set_xlim([xmin,xmax])
    ax0.set_ylim([0,dosmax*1.1])
    ax0.tick_params(axis="x", bottom=True, top=True, direction="in", color=colpal[2], labelcolor=colpal[4], width=1, zorder=0)
else : # vertical
    fig=plt.figure(figsize=(1,3.75))
    if doscar0.Ns==1 :
        gs0=fig.add_gridspec(1,1,wspace=0.0,hspace=0.00,left=0.02,right=0.98,top=0.97, bottom=0.07)
        ax0=gs0.subplots()
    else :
        gs0=fig.add_gridspec(1,2,wspace=0.0,hspace=0.00,left=0.02,right=0.98,top=0.97, bottom=0.07)
        (ax1,ax0)=gs0.subplots()

    ax0.axhline(linewidth=1,color=colpal[2],zorder=0)
    ax0.plot(doscar0.dos[0],doscar0.energy,color=colpal[0],linewidth=1,zorder=3)
    ax0.set_xlim([0,dosmax*1.1])
    ax0.set_ylim([xmin,xmax])
    ax0.tick_params(axis="x", bottom=False, top=False, direction="in", length=0)
    
ax0.tick_params(axis="y", left=True, right=True, direction="in", color=colpal[2], labelcolor=colpal[4], width=1, zorder=0)
for edge in ["bottom","top","left","right"] :
    ax0.spines[edge].set_color(colpal[4])
    ax0.spines[edge].set_linewidth(1)
    ax0.spines[edge].set_zorder(4)

if fvertical==False:
    if doscar0.Ns==1 :
        ax0.set_xlabel("Energy (eV)",color=colpal[4],labelpad=-1)
        ax0.set_ylabel("DOS (eV⁻¹)",color=colpal[4])
    
    elif doscar0.Ns==2 : 
        ax0.set_xticklabels([])
        ax1.axvline(linewidth=1,color=colpal[2],zorder=0)
        ax1.plot(doscar0.energy,doscar0.dos[1],color=colpal[1],linewidth=1,zorder=3)
        ax1.set_xlim([xmin,xmax])
        ax1.set_ylim([dosmax*1.1,0.0])
        ax1.set_xlabel("Energy (eV)",color=colpal[4],labelpad=-1)
        ax0.set_ylabel("DOS (eV⁻¹)",color=colpal[4],y=0.0)
        ax1.tick_params(axis="x", bottom=True, top=True, direction="in", color=colpal[2], labelcolor=colpal[4], width=1, zorder=0)
else: # vertical
    ax0.set_yticklabels([])
    if doscar0.Ns==1 :
        # put "DOS" as ticklabel to align with the K point labels in bs
        ax0.set_xticks([dosmax*0.55],["DOS"],color=colpal[4])
    
    elif doscar0.Ns==2 : 
        ax0.set_xticks([0.0],["DOS"],color=colpal[4])
        ax1.set_xticks([],[])
        ax1.set_yticklabels([])
        ax1.axhline(linewidth=1,color=colpal[2],zorder=0)
        ax1.plot(doscar0.dos[1],doscar0.energy,color=colpal[1],linewidth=1,zorder=3)
        ax1.set_xlim([dosmax*1.1,0.0])
        ax1.set_ylim([xmin,xmax])
        ax1.tick_params(axis="x", bottom=False, top=False, direction="in", length=0)

if doscar0.Ns==2: 
    ax1.tick_params(axis="y", left=True, right=True, direction="in", color=colpal[2], labelcolor=colpal[4], width=1, zorder=0)
    for edge in ["bottom","top","left","right"] :
        ax1.spines[edge].set_color(colpal[4])
        ax1.spines[edge].set_linewidth(1)
        ax1.spines[edge].set_zorder(4)

fig.savefig(filename_out,dpi=1200)
