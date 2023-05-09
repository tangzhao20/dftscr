#!/usr/bin/env python3

# bands.py package (E1) (E2)
# Making the band structure plot

# if E2 exists, the energy range is (E1,E2)
# if only E1 exists, the energy range is (-E1,E1)
# if neither exists, the energy range is (-5 eV, 5 eV)

import os
import sys
from classes import poscar, eigenval, doscar, kpoints_band, procar
from commons import load_packagename, load_palette
import matplotlib.pyplot as plt
import matplotlib as mpl

fsecond=False

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
        eigenval1.eigshift(eigenval1.vbm)
        eigenval1.writegap(kpoints1)
    else :
        doscar1=doscar.DOSCAR()
        eigenval1.eigshift(doscar1.ef)
    
    x=eigenval1.bandkpout(kp=kpoints1,reclc=rlc)
    eigout=eigenval1.eigtrans()
    
    kphx=kpoints1.kphx_out(rlc)
    kphlabel=kpoints1.kphlabel_out()


elif package in packagename["qe"]+packagename["qeproj"] :
    # Input: *.xml, nscf.in, kpath.in
    # No need to run bands.x

    poscar1=poscar.POSCAR(empty=True)
    poscar1.fileread_qe("nscf.in")
    rlc=poscar1.reclc_out()
    
    eigenval1=eigenval.EIGENVAL(empty=True)
    eigenval1.fileread_qexml()
    
    kpoints1=kpoints_band.KPOINTS_band(empty=True)
    kpoints1.fileread_kpathin()
    
    if eigenval1.is_semic==True :
        eigenval1.eigshift(eigenval1.vbm)
        eigenval1.writegap(kpoints1)
    else :
        doscar1=doscar.DOSCAR(empty=True)
        doscar1.fileread_xml(filename)
        eigenval1.eigshift(doscar1.ef)
    
    x=eigenval1.bandkpout(kp=kpoints1,reclc=rlc)
    print("x=",x)
    eigout=eigenval1.eigtrans()
    
    kphx=kpoints1.kphx_out(rlc)
    kphlabel=kpoints1.kphlabel_out()
    
    ## the x axis somehow doesn`t match here. qe used the unit 2pi/alat?
    #factor=max(kphx)/max(x)
    #for ix in range(len(x)): 
    #    x[ix]=x[ix]*factor

elif package in packagename["wannier90"] :
    # Input : nscf.in, ../bands/*.xml, *_band.kpt, *_band.dat, kpath.in
    # Only semiconductors are supported

    poscar1=poscar.POSCAR(empty=True)
    poscar1.fileread_qe("nscf.in")
    rlc=poscar1.reclc_out()

    pad=0
    for w in sys.argv :
        if w.startswith("pad=") :
            pad=int(w.split("=")[1])
            sys.argv.remove(w)
            break

    eigenval1=eigenval.EIGENVAL(empty=True)
    eigenval1.fileread_wan(Nb_pad=pad)

    kpoints1=kpoints_band.KPOINTS_band(empty=True)
    kpoints1.fileread_kpathin()

    # find a ../bands/*.xml file
    files = os.listdir("../bands")
    for f in files:
        if f.endswith('.xml'):
            filename=f
            fsecond=True
            break

    if fsecond :
        eigenval2=eigenval.EIGENVAL(empty=True)
        eigenval2.fileread_qexml("../bands/"+filename)

        eigenval2.gap()
        eigenval2.eigshift(eigenval2.vbm)

    eigenval1.occ=[]
    for ik in range(eigenval1.Nk) :
        occ0=[]
        for ib in range(eigenval2.Nvb[0]-pad) :
            occ0.append([1.0])
        if fsecond :
           for ib in range(eigenval2.Nvb[0]-pad,eigenval1.Nb) :
               occ0.append([0.0])
           eigenval1.occ.append(occ0)

    eigenval1.gap()
    eigenval1.eigshift(eigenval1.vbm)

    x=eigenval1.bandkpout(kp=kpoints1,reclc=rlc)
    eigout=eigenval1.eigtrans()

    kphx=kpoints1.kphx_out(rlc)
    kphlabel=kpoints1.kphlabel_out()

    if fsecond :
        x2=eigenval2.bandkpout(kp=kpoints1)
        eigout2=eigenval2.eigtrans()

        # the x axis somehow doesn`t match here. qe used the unit 2pi/alat?
        factor=max(kph)/max(x2)
        for ix in range(len(x2)) :
            x2[ix]=x2[ix]*factor

else:
    print("Package \""+package+"\" is not supported yet.")
    print("python bands.py package (Emin) (Emax)")
    sys.exit()

lproj=False
if package in packagename["vaspproj"]+packagename["qeproj"] :
    lproj=True
    if eigenval1.Ns==2 :
        print("Spin polarized projected band structure not supported yet.\n")
        sys.exit()
    if package in packagename["vaspproj"] :
        # Input: PROCAR
        procar1=procar.PROCAR()
    elif package in packagename["qeproj"] :
        # Input: projwfc.out
        procar1=procar.PROCAR(empty=True)
        procar1.fileread_qe()
    
    atomlist=sys.argv[2]
    del sys.argv[2]
    atomflag=procar1.readatomlist(atomlist,poscar1)
    orblist=sys.argv[2]
    del sys.argv[2]
    orbflag=procar1.readorblist(orblist)

if len(sys.argv)>=4 :
    ymax=float(sys.argv[3])
    ymin=float(sys.argv[2])
elif len(sys.argv)==3 :
    ymax=float(sys.argv[2])
    ymin=-float(sys.argv[2])
else :
    ymax=5.0
    ymin=-5.0

colpal=load_palette() #color palette: [blue, orange, gray, white, black]
mpl.rcParams["font.sans-serif"].insert(0,"Noto Sans")
mpl.rcParams.update({'font.size': 14})

# band structure plot

fig=plt.figure(figsize=(5,3.75))
gs0=fig.add_gridspec(1,1,wspace=0.0,hspace=0.00,left=0.14,right=0.98,top=0.97, bottom=0.07)
ax0=gs0.subplots()

ax0.grid(axis="x",linewidth=1, color=colpal[2],zorder=0)
for il in range(len(kphlabel)):
    if "|" in kphlabel[il]:
        ax0.axvline(x=kphx[il],linewidth=1,color=colpal[4],zorder=4)
ax0.axhline(linewidth=1,color=colpal[2],zorder=0)

spinlabel=["spin up","spin down"]
for s in range(eigenval1.Ns) :
    for b in range(eigenval1.Nb) :
        if eigenval1.Ns==2 and b==0 :
            ax0.plot(x,eigout[s][b],color=colpal[s],label=spinlabel[s],linewidth=1,zorder=3-s)
        else :
            ax0.plot(x,eigout[s][b],color=colpal[s],linewidth=1,zorder=3-s)

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
    outputname="proj_"+atomlist+"_"+orblist+".png"
    
    f2=open("proj_"+atomlist+"_"+orblist+".dat","w")
    f2.write("#k-point energy(eV) pointsize\n")
    for ib in range(eigenval1.Nb):
        for ik in range(eigenval1.Nk):
            f2.write(str(x[ik])+" "+str(eigout[0][ib][ik])+" "+str(projplotsize[ib][ik])+"\n")
        f2.write("\n")
    f2.close()

# second band structure plot (for wannier)

if fsecond :
    for b in range(eigenval2.Nb) :
        ax0.plot(x2,eigout2[0][b],color=colpal[1],linewidth=1,zorder=2)
    f3=open("eigenval2.dat","w")
    for i in range(len(eigout2[0])):
        for j in range(len(x2)) :
            f3.write(str(x2[j])+" "+str(eigout2[0][i][j])+" "+"\n")
        f3.write("\n")
    f3.close()

ax0.set_ylim(ymin,ymax)
ax0.set_xlim(x[0],x[len(x)-1])
ax0.set_xticks(kphx,kphlabel,color=colpal[4])
ax0.tick_params(axis="x", direction="in", length=0)
ax0.tick_params(axis="y", left=True, right=True, direction="in", color=colpal[2], labelcolor=colpal[4], width=1, zorder=0)
ax0.set_ylabel("Energy (eV)",labelpad=-2,color=colpal[4])
for edge in ["bottom", "top", "left", "right"] :
    ax0.spines[edge].set_color(colpal[4])
    ax0.spines[edge].set_linewidth(1)
    ax0.spines[edge].set_zorder(4)

f4=open("label.dat","w")
for i in range(len(kphx)):
    f4.write(str(kphx[i])+" "+kphlabel[i]+"\n")
f4.close()

fig.savefig(outputname,dpi=1200)

