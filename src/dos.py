#!/usr/bin/env python3

# print DOS

# Input: DOSCAR

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from classes import doscar
from commons import load_packagename, load_palette

if len(sys.argv)>=3 :
    xmax=float(sys.argv[2])
    xmin=float(sys.argv[1])
elif len(sys.argv)==2 :
    xmax=float(sys.argv[1])
    xmin=-float(sys.argv[1])
else :
    xmax=5.0
    xmin=-5.0

doscar0=doscar.DOSCAR()

colpal=load_palette() #color palette: [blue, orange, gray, white, black]
mpl.rcParams["font.sans-serif"].insert(0,"Noto Sans")
mpl.rcParams.update({'font.size': 14})

fig=plt.figure(figsize=(5,3.75))
if doscar0.Ns==1 :
    gs0=fig.add_gridspec(1,1,wspace=0.0,hspace=0.00,left=0.17,right=0.95,top=0.97, bottom=0.15)
    ax0=gs0.subplots()
else :
    gs0=fig.add_gridspec(2,1,wspace=0.0,hspace=0.00,left=0.17,right=0.95,top=0.97, bottom=0.15)
    (ax0,ax1)=gs0.subplots()

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
        
ax0.axvline(linewidth=1,color=colpal[2],zorder=0)
ax0.plot(doscar0.energy,doscar0.dos[0],color=colpal[0],linewidth=1,zorder=3)
ax0.set_xlim([xmin,xmax])
ax0.set_ylim([0,dosmax*1.1])
ax0.tick_params(axis="x", bottom=True, top=True, direction="in", color=colpal[2], labelcolor=colpal[4], width=1, zorder=0)
ax0.tick_params(axis="y", direction="in", color=colpal[2], labelcolor=colpal[4], width=1, zorder=0)

for edge in ["bottom","top","left","right"] :
    ax0.spines[edge].set_color(colpal[4])
    ax0.spines[edge].set_linewidth(1)
    ax0.spines[edge].set_zorder(4)

if doscar0.Ns==1 :
    ax0.set_xlabel("Energy (eV)",color=colpal[4])
    ax0.set_ylabel("DOS (eV⁻¹)",color=colpal[4])

elif doscar0.Ns==2 : 
    ax0.set_xticklabels([])
    ax1.axvline(linewidth=1,color=colpal[2],zorder=0)
    ax1.plot(doscar0.energy,doscar0.dos[1],color=colpal[1],linewidth=1,zorder=3)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([dosmax*1.1,0.0])
    ax1.tick_params(axis="x", bottom=True, top=True, direction="in", color=colpal[2], labelcolor=colpal[4], width=1, zorder=0)
    ax1.tick_params(axis="y", direction="in", color=colpal[2], labelcolor=colpal[4], width=1, zorder=0)
    ax1.set_xlabel("Energy (eV)",color=colpal[4])
    ax0.set_ylabel("DOS (eV⁻¹)",color=colpal[4],y=0.0)
    for edge in ["bottom","top","left","right"] :
        ax1.spines[edge].set_color(colpal[4])
        ax1.spines[edge].set_linewidth(1)
        ax1.spines[edge].set_zorder(4)

fig.savefig("dos.png",dpi=1200)
