#!/bin/python3
import sys
#sys.path.append("/user/ztang5/codes/headers")
import dftscr.vaspfiles as vaspfiles
#import vaspfiles
import matplotlib.pyplot as plt

poscar1=vaspfiles.poscar.POSCAR()
#poscar1=filein.POSCAR("CONTCAR")
rlc=poscar1.reclc_out()

eigenval1=vaspfiles.eigenval.EIGENVAL()
kpoints1=vaspfiles.kpoints_band.KPOINTS_band()

if eigenval1.is_semic==True :
    vbm=eigenval1.vbm
    eigenval1.eigshift(vbm)
    eigenval1.writegap(kpoints1)
else :
    doscar1=vaspfiles.doscar.DOSCAR()
    ef=doscar1.ef_out()
    eigenval1.eigshift(ef)

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
plt.ylim(ymin,ymax)
plt.xlim(x[0],x[len(x)-1])
plt.xticks(kphout,kphlabel)

f2=open("label.dat","w")
for i in range(len(kphout)):
    f2.write(str(kphout[i])+" "+kphlabel[i]+"\n")
f2.close()
plt.ylabel("Energy (eV)")
plt.savefig("bs.png",dpi=300)
#plt.show()
           
 
