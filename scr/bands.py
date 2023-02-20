#!/bin/python3

# bands.py package (E1) (E2)
# Making the band structure plot

# if E2 exists, the energy range is (E1,E2)
# if only E1 exists, the energy range is (-E1,E1)
# if neither exists, the energy range is (-5 eV, 5 eV)

import os
import sys
from dftscr.vaspfiles import poscar, eigenval, doscar, kpoints_band, procar
from commons import load_packagename
import matplotlib.pyplot as plt
import matplotlib as mpl

if len(sys.argv)<=1 :
    print("python bands.py package (Emin) (Emax)")
    sys.exit()
package=sys.argv[1]

packagename=load_packagename()

if package in packagename["vasp"]+packagename["vaspproj"] : 
    # Input: EIGENVAL, KPOINTS, POSCAR, (DOSCAR)

    poscar1=poscar.POSCAR()
    rlc=poscar1.reclc_out()
    
    eigenval1=eigenval.EIGENVAL()
    kpoints1=kpoints_band.KPOINTS_band()
    
    if eigenval1.is_semic==True :
        vbm=eigenval1.vbm
        eigenval1.eigshift(vbm)
        eigenval1.writegap(kpoints1)
    else :
        doscar1=doscar.DOSCAR()
        ef=doscar1.ef
        eigenval1.eigshift(ef)
    
    x=eigenval1.bandkpout(reclc=rlc)
    eigout=eigenval1.eigtrans()
    
    kphout=kpoints1.kphout_out(rlc)
    kphlabel=kpoints1.kph_out()


elif package in packagename["qe"] :
    # Input: *.xml, nscf.in, kpath.in

    poscar1=poscar.POSCAR(empty=True)
    poscar1.fileread_qe("nscf.in")
    rlc=poscar1.reclc_out()
    
    eigenval1=eigenval.EIGENVAL(empty=True)
    # find a .xml file
    files = os.listdir()
    for f in files:
        if f.endswith('.xml'):
            filename=f
            break
    eigenval1.fileread_qexml(filename)
    
    kpoints1=kpoints_band.KPOINTS_band(empty=True)
    kpoints1.fileread_kpathin()
    
    if eigenval1.is_semic==True :
        vbm=eigenval1.vbm
        eigenval1.eigshift(vbm)
        eigenval1.writegap(kpoints1)
    else :
        doscar1=doscar.DOSCAR(empty=True)
        doscar1.fileread_xml(filename)
        ef=doscar1.ef
        eigenval1.eigshift(ef)
    
    x=eigenval1.bandkpout()
    eigout=eigenval1.eigtrans()
    
    kphout=kpoints1.kphout_out(rlc)
    kphlabel=kpoints1.kph_out()
    
    # the x axis somehow doesn`t match here. qe used the unit 2pi/alat?
    factor=max(kphout)/max(x)
    for ix in range(len(x)): 
        x[ix]=x[ix]*factor

else:
    print("Package \""+package+"\" is not supported yet.")
    print("python bands.py package (Emin) (Emax)")
    sys.exit()

lproj=False
if package in packagename["vaspproj"] :
    # Input: PROCAR

    lproj=True
    if eigenval1.Ns==2 :
        print("Spin polarized projected band structure not supported yet.\n")
        sys.exit()
    procar1=procar.PROCAR()
    atomlist=sys.argv[2]
    orblist=sys.argv[3]
    orbflag=procar1.readorblist(orblist)
    atomflag=procar1.readatomlist(atomlist,poscar1)

if package in packagename["vasp"]+packagename["qe"] :
    if len(sys.argv)>=4 :
        ymax=float(sys.argv[3])
        ymin=float(sys.argv[2])
    elif len(sys.argv)==3 :
        ymax=float(sys.argv[2])
        ymin=-float(sys.argv[2])
    else :
        ymax=5.0
        ymin=-5.0
elif package in packagename["vaspproj"] :
    if len(sys.argv)>=6 :
        ymax=float(sys.argv[5])
        ymin=float(sys.argv[4])
    elif len(sys.argv)==5 :
        ymax=float(sys.argv[4])
        ymin=-float(sys.argv[4])
    else :
        ymax=5.0
        ymin=-5.0
else :
    print("Package \""+package+"\" is not supported yet.")
    print("This occurs when reading the y range. Check code.")
    print("python bands.py package (Emin) (Emax)")
    sys.exit()

colpal=["#005f86","#f8971f","#d6d2c4","#ffffff","#333f48"] #color palette: [Blue, Orange, Gray, white, black]
mpl.rcParams["font.sans-serif"].insert(0,"Noto Sans")

# band structure plot

fig=plt.figure(figsize=(4,3))
gs0=fig.add_gridspec(1,1,wspace=0.0,hspace=0.00,left=0.15,right=0.98,top=0.97, bottom=0.07)
ax0=gs0.subplots()

ax0.grid(axis="x",linewidth=1, color=colpal[2],zorder=0)
ax0.axhline(linewidth=1,color=colpal[2],zorder=0)

spinlabel=["spin up","spin down"]
for s in range(eigenval1.Ns) :
    for b in range(eigenval1.Nb) :
        if eigenval1.Ns==2 and b==0 :
            ax0.plot(x,eigout[s][b],color=colpal[s],label=spinlabel[s],linewidth=1,zorder=3)
        else :
            ax0.plot(x,eigout[s][b],color=colpal[s],linewidth=1,zorder=3)

outputname="bs.png"

#if eigenval1.Ns==2 :
#    plt.legend()

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

# projection plot

if lproj:
    dotsize=50.0
    projplotsize=procar1.plot(atomflag,orbflag,dotsize)
    for ib in range(eigenval1.Nb):
        ax0.scatter(x,eigout[0][ib],s=projplotsize[ib],c=colpal[1],zorder=2)
    outputname="proj_"+sys.argv[2]+"_"+sys.argv[3]+".png"
    
    f2=open("proj_"+sys.argv[2]+"_"+sys.argv[3]+".dat","w")
    f2.write("#k-point energy(eV) pointsize\n")
    for ib in range(eigenval1.Nb):
        for ik in range(eigenval1.Nk):
            f2.write(str(x[ik])+" "+str(eigout[0][ib][ik])+" "+str(projplotsize[ib][ik])+"\n")
        f2.write("\n")
    f2.close()

ax0.set_ylim(ymin,ymax)
ax0.set_xlim(x[0],x[len(x)-1])
ax0.set_xticks(kphout,kphlabel,color=colpal[4])
ax0.tick_params(axis="x", length=0)
ax0.tick_params(axis="y", colors=colpal[4])
ax0.tick_params(axis="y", direction="in", color=colpal[2], width=1,zorder=0)
ax0.set_ylabel("Energy (eV)",labelpad=-2,color=colpal[4])
ax0.spines['bottom'].set_color(colpal[4])
ax0.spines['top'].set_color(colpal[4])
ax0.spines['left'].set_color(colpal[4])
ax0.spines['right'].set_color(colpal[4])
ax0.spines['bottom'].set_linewidth(1)
ax0.spines['top'].set_linewidth(1)
ax0.spines['left'].set_linewidth(1)
ax0.spines['right'].set_linewidth(1)
ax0.spines['bottom'].set_zorder(4)
ax0.spines['top'].set_zorder(4)
ax0.spines['left'].set_zorder(4)
ax0.spines['right'].set_zorder(4)

f4=open("label.dat","w")
for i in range(len(kphout)):
    f4.write(str(kphout[i])+" "+kphlabel[i]+"\n")
f4.close()

fig.savefig(outputname,dpi=1200)

