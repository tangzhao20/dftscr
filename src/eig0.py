#!/usr/bin/env python3

# eig0.py package (E1) (E2)
# Plot the eigenvalues at the Gamma point.

# if E2 exists, the energy range is (E1,E2)
# if only E1 exists, the energy range is (-E1,E1)
# if neither exists, the energy range is (-5 eV, 5 eV)

import os
import sys
from classes import EIGENVAL
from commons import load_packagename, load_palette
import matplotlib.pyplot as plt
import matplotlib as mpl

fsecond=False

if len(sys.argv)<=1 :
    print("python bands.py package (Emin) (Emax)")
    sys.exit()
package=sys.argv[1]

packagename=load_packagename()

eigenval1=EIGENVAL()
if package in packagename["vasp"] : 
    # Input: EIGENVAL
    # TODO: test VASP
    eigenval1.fileread_vasp()
    
elif package in packagename["qe"] :
    # Input: *.xml
    eigenval1.fileread_qexml()

elif package in packagename["parsec"] :
    # Input: bands.dat
    eigenval1.fileread_parsec()

else:
    print("Package \""+package+"\" is not supported yet.")
    print("python bands.py package (Emin) (Emax)")
    sys.exit()

if eigenval1.is_semic==True :
    eigenval1.eigshift(eigenval1.vbm)
else :
    print("Error: only semiconductor with bandgap is supported")
    sys.exit()

ik_gamma=-1
for ik in range(eigenval1.Nk) :
    summ=0
    for ix in range(3) :
        summ+=abs(eigenval1.kp[ik][ix])
    if summ<1e-6 :
        ik_gamma=ik
        break
if ik_gamma==-1 :
    print("Error: The Gamma point is not found")
    sys.exit()

if len(sys.argv)>=4 :
    ymax=float(sys.argv[3])
    ymin=float(sys.argv[2])
elif len(sys.argv)==3 :
    ymax=float(sys.argv[2])
    ymin=-float(sys.argv[2])
else :
    ymax=5.0
    ymin=-5.0

palette=load_palette()
mpl.rcParams["font.sans-serif"].insert(0,"Noto Sans")
mpl.rcParams.update({'font.size': 14})

# band structure plot

linecolor=["darkblue","orange"]

fig=plt.figure(figsize=(1.8,3.75))
gs0=fig.add_gridspec(1,eigenval1.Ns,wspace=0.0,hspace=0.00,left=0.45,right=0.97,top=0.97, bottom=0.07)

ax=[]
for ispin in range(eigenval1.Ns):
    ax.append(fig.add_subplot(gs0[ispin]))

    ax[ispin].grid(axis="x",linewidth=1, color=palette["gray"],zorder=0)
    ax[ispin].axhline(linewidth=1,color=palette["gray"],zorder=0)

    for ib in range(eigenval1.Nb) :
        if eigenval1.occ[ik_gamma][ib][ispin] >0.5 :
            linestyle="solid"
        else :
            linestyle="dashed"
        ax[ispin].axhline(y=eigenval1.eig[ik_gamma][ib][ispin],linestyle=linestyle,linewidth=1,color=palette[linecolor[ispin]],zorder=2)

outputname="eig0.png"

ax[0].set_ylabel("Energy (eV)",labelpad=-2,color=palette["black"])
for ispin in range(eigenval1.Ns):
    ax[ispin].set_ylim(ymin,ymax)
    ax[ispin].set_xlim(0,1)
    ax[ispin].set_xticks([])
    ax[ispin].tick_params(axis="x", direction="in", length=0)
    ax[ispin].tick_params(axis="y", left=False, right=False, direction="in", color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
    if ispin!=0 :
        ax[ispin].yaxis.set_ticklabels([])
    for edge in ["bottom", "top", "left", "right"] :
        ax[ispin].spines[edge].set_color(palette["black"])
        ax[ispin].spines[edge].set_linewidth(1)
        ax[ispin].spines[edge].set_zorder(4)
ax[0].tick_params(axis="y", left=True)
ax[-1].tick_params(axis="y", right=True)

fig.savefig(outputname,dpi=1200)

