#!/bin/python3

# eigenval.py (E1) (E2)
# Making the band structure plot

# if E2 exists, the energy range is (E1,E2)
# if only E1 exists, the energy range is (-E1,E1)
# if neither exists, the energy range is (-5 eV, 5 eV)

# Input: EIGENVAL, KPOINTS, POSCAR (DOSCAR)

import sys
from dftscr.vaspfiles import poscar, eigenval, doscar, kpoints_band
import matplotlib.pyplot as plt

poscar1=poscar.POSCAR()
#poscar1=poscar.POSCAR("CONTCAR")
rlc=poscar1.reclc_out()

eigenval1=eigenval.EIGENVAL()
kpoints1=kpoints_band.KPOINTS_band()

if eigenval1.is_semic==True :
    vbm=eigenval1.vbm
    eigenval1.eigshift(vbm)
    eigenval1.writegap(kpoints1)
else :
    doscar1=doscar.DOSCAR()
    ef=doscar1.ef_out()
    eigenval1.eigshift(ef)

x=eigenval1.bandkpout(rlc)
eigout=eigenval1.eigtrans()

kphout=kpoints1.kphout_out(rlc)
kphlabel=kpoints1.kphl_out()

# make plot

fig=plt.figure(figsize=(4,3))
gs0=fig.add_gridspec(1,1,wspace=0.0,hspace=0.00,left=0.15,right=0.95,top=0.95, bottom=0.15)
ax0=gs0.subplots()

ax0.grid(axis="x",linewidth=0.75, color="silver")
ax0.axhline(linewidth=0.75,color="silver")

spinlabel=["spin up","spin down"]
for s in range(eigenval1.Ns) :
    for b in range(eigenval1.Nb) :
        if eigenval1.Ns==2 and b==0 :
            ax0.plot(x,eigout[s][b],color="C"+str(s),label=spinlabel[s])
        else :
            ax0.plot(x,eigout[s][b],color="C"+str(s))

f3=open("eigenval.dat","w")
if eigenval1.Ns==2 :
    for i in range(len(eigout[0])):
        for j in range(len(x)) :
            f3.write(str(x[j])+" "+str(eigout[0][i][j])+" "+str(eigout[1][i][j])+"\n")
        f3.write("\n")
else :
    for i in range(len(eigout[0])):
        for j in range(len(x)) :
            f3.write(str(x[j])+" "+str(eigout[0][i][j])+" "+"\n")
        f3.write("\n")
f3.close()
#if eigenval1.Ns==2 :
#    plt.legend()

if len(sys.argv)>=3 :
    ymax=float(sys.argv[2])
    ymin=float(sys.argv[1])
elif len(sys.argv)==2 :
    ymax=float(sys.argv[1])
    ymin=-float(sys.argv[1])
else :
    ymax=5.0
    ymin=-5.0
ax0.set_ylim(ymin,ymax)
ax0.set_xlim(x[0],x[len(x)-1])
ax0.set_xticks(kphout,kphlabel)
ax0.set_ylabel("Energy (eV)")

f2=open("label.dat","w")
for i in range(len(kphout)):
    f2.write(str(kphout[i])+" "+kphlabel[i]+"\n")
f2.close()

fig.savefig("bs.png",dpi=300)
#plt.show()
           
 
