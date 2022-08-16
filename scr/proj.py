#!/bin/python3

# python proj.py atoms orbitals

import sys
#sys.path.append("/user/ztang5/codes/headers")
from dftscr.vaspfiles import *
import matplotlib.pyplot as plt

poscar1=poscar.POSCAR("CONTCAR")
rlc=poscar1.reclc_out()

doscar1=doscar.DOSCAR()
ef=doscar1.ef_out()

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

if eigenval1.Ns==2 :
    print("Spin polarized projected band structure not supported yet.\n")
    sys.exit()

procar1=procar.PROCAR()
atomlist=sys.argv[1]
orblist=sys.argv[2]
orbflag=procar1.readorblist(orblist)
atomflag=procar1.readatomlist(atomlist,poscar1)

x=eigenval1.bandkpout(rlc)
eigout=eigenval1.eigtrans()

kphout=kpoints1.kphout_out(rlc)
kphlabel=kpoints1.kphl_out()

# make plot
plt.gca().grid(axis="x",linewidth=0.75, color="silver")
plt.axhline(linewidth=0.75,color="silver")

spinlabel=["spin up","spin down"]
for s in range(eigenval1.Ns) :
    for b in range(eigenval1.Nb) :
        if eigenval1.Ns==2 and b==0 :
            plt.plot(x,eigout[s][b],color="C"+str(s),label=spinlabel[s])
        else :
            plt.plot(x,eigout[s][b],color="C"+str(s))


fac=50.0
projplotsize=procar1.plot(atomflag,orbflag,fac)
for ib in range(eigenval1.Nb):
    plt.scatter(x,eigout[0][ib],s=projplotsize[ib],c="C1")

f2=open("proj_"+sys.argv[1]+"_"+sys.argv[2]+".dat","w")
f2.write("#k-point energy(eV) pointsize\n")
for ib in range(eigenval1.Nb):
    for ik in range(eigenval1.Nk):
        f2.write(str(x[ik])+" "+str(eigout[0][ib][ik])+" "+str(projplotsize[ib][ik])+"\n")
    f2.write("\n")    
f2.close()
    
if eigenval1.Ns==2 :
    plt.legend()

if len(sys.argv)>=5:
    ymax=float(sys.argv[4])
    ymin=float(sys.argv[3])
elif len(sys.argv)==2 :
    ymax=float(sys.argv[3])
    ymin=-float(sys.argv[3])
else :
    ymax=5.0
    ymin=-5.0
plt.ylim(ymin,ymax)
plt.xlim(x[0],x[len(x)-1])

plt.xticks(kphout,kphlabel)
f3=open("xlabel.dat","w")
for ik in range(len(kphout)):
    f3.write(str(kphout[ik])+" "+kphlabel[ik]+"\n")
f3.close()

plt.ylabel("Energy (eV)")
plt.savefig("proj_"+sys.argv[1]+"_"+sys.argv[2]+".png",dpi=300,bbox_inches='tight')
#plt.show()
           
 
