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
    # TODO: test VASP

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
    energy=eigenval1.eigtrans()
    
    xticks=kpoints1.xticks_out(rlc)
    xlabels=kpoints1.xlabels_out()

elif package in packagename["qe"]+packagename["qeproj"] :
    # Input: *.xml, kpath.in
    # No need to run bands.x

    poscar1=poscar.POSCAR(empty=True)
    #poscar1.fileread_qe("nscf.in")
    poscar1.fileread_xml()
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
        doscar1.fileread_xml()
        eigenval1.eigshift(doscar1.ef)
    
    x=eigenval1.bandkpout(kp=kpoints1,reclc=rlc)
    energy=eigenval1.eigtrans()
    
    xticks=kpoints1.xticks_out(rlc)
    xlabels=kpoints1.xlabels_out()
    
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
        if eigenval2.is_semic==True :
            eigenval1.is_semic=True

    eigenval1.occ=[]
    for ik in range(eigenval1.Nk) :
        occ0=[]
        for ib in range(eigenval2.Nvb[0]-pad) :
            occ0.append([1.0])
        if fsecond :
           for ib in range(eigenval2.Nvb[0]-pad,eigenval1.Nb) :
               occ0.append([0.0])
           eigenval1.occ.append(occ0)

    if eigenval1.is_semic==True :
        eigenval1.gap()
        eigenval1.writegap(kpoints1)
        eigenval1.eigshift(eigenval1.vbm)
    else : 
        print("Metal band structure are not shifted")

    x=eigenval1.bandkpout(kp=kpoints1,reclc=rlc)
    energy=eigenval1.eigtrans()

    xticks=kpoints1.xticks_out(rlc)
    xlabels=kpoints1.xlabels_out()

    if fsecond :
        x2=eigenval2.bandkpout(kp=kpoints1,reclc=rlc)
        energy2=eigenval2.eigtrans()

else:
    print("Package \""+package+"\" is not supported yet.")
    print("python bands.py package (Emin) (Emax)")
    sys.exit()

# width is the physical width of each plots
# Nx is the number of x points in each plots
# ixl and ixr are the left and right index of each panel, to seperate the bands
width=[]
Nx=[]
ixl=[]
ixr=[]
for p in x :
    width.append(max(p))
    Nx.append(len(p))
    if ixl==[] :
        ixl=[0]
    else :
        ixl.append(ixr[-1])
    ixr.append(ixl[-1]+Nx[-1])
if fsecond :
    Nx2=[]
    ix2l=[]
    ix2r=[]
    for p in x2 :
        Nx2.append(len(p))
        if ix2l==[] :
            ix2l=[0]
        else :
            ix2l.append(ix2r[-1])
        ix2r.append(ix2l[-1]+Nx2[-1])

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
gs0=fig.add_gridspec(1,len(xticks),wspace=0.0,hspace=0.00,left=0.14,right=0.98,top=0.97, bottom=0.07,width_ratios=width[:len(xticks)])
ax=[]
for ip in range(len(xticks)):
    ax.append(fig.add_subplot(gs0[ip]))

    ax[ip].grid(axis="x",linewidth=1, color=colpal[2],zorder=0)
    ax[ip].axhline(linewidth=1,color=colpal[2],zorder=0)
    spinlabel=["spin up","spin down"]
    for s in range(eigenval1.Ns) :
        for b in range(eigenval1.Nb) :
            if eigenval1.Ns==2 and b==0 :
                ax[ip].plot(x[ip],energy[s][b][ixl[ip]:ixr[ip]],color=colpal[s],label=spinlabel[s],linewidth=1,zorder=3-s)
            else :
                ax[ip].plot(x[ip],energy[s][b][ixl[ip]:ixr[ip]],color=colpal[s],linewidth=1,zorder=3-s)

outputname="bs.png"

#if eigenval1.Ns==2 :
#    plt.legend()

f3=open("eigenval.dat","w")
if eigenval1.Ns==2 :
    for ib in range(len(energy[0])):
        for ip in range(len(xticks)):
            ik0=0
            for ik in range(len(x[ip])):
                f3.write(str(x[ip][ik])+" "+str(energy[0][ib][ik0])+" "+str(energy[1][ib][ik0])+"\n")
                ik0+=1
        f3.write("\n")
else :
    for ib in range(len(energy[0])):
        for ip in range(len(xticks)):
            ik0=0
            for ik in range(len(x[ip])):
                f3.write(str(x[ip][ik])+" "+str(energy[0][ib][ik0])+" "+"\n")
                ik0+=1
        f3.write("\n")
f3.close()

# projection plot

if lproj:
    dotsize=50.0
    projplotsize=procar1.plot(atomflag,orbflag,dotsize)
    for ip in range(len(xticks)):
        for ib in range(eigenval1.Nb):
            ax[ip].scatter(x[ip],energy[0][ib][ixl[ip]:ixr[ip]],s=projplotsize[ib][ixl[ip]:ixr[ip]],c=colpal[1],zorder=2)
    outputname="proj_"+atomlist+"_"+orblist+".png"
    
    f2=open("proj_"+atomlist+"_"+orblist+".dat","w")
    f2.write("#k-point energy(eV) pointsize\n")
    for ib in range(eigenval1.Nb):
        for ip in range(len(x)):
            ik0=0
            for ik in range(len(x[ip])):
                f2.write(str(x[ip][ik])+" "+str(energy[0][ib][ik0])+" "+str(projplotsize[ib][ik0])+"\n")
                ik0+=1
        f2.write("\n")
    f2.close()

# second band structure plot (for wannier)

if fsecond :
    for ip in range(len(xticks)):
        for b in range(eigenval2.Nb):
            ax[ip].plot(x2[ip],energy2[0][b][ix2l[ip]:ix2r[ip]],color=colpal[1],linewidth=1,zorder=2)
    f3=open("eigenval2.dat","w")
    for ib in range(len(energy2[0])):
        for ip in range(len(xticks)):
            ik0=0
            for ik in range(len(x2[ip])) :
                f3.write(str(x2[ip][ik])+" "+str(energy2[0][ib][ik0])+" "+"\n")
                ik0+=1
        f3.write("\n")
    f3.close()

ax[0].set_ylabel("Energy (eV)",labelpad=-2,color=colpal[4])
for ip in range(len(xticks)):
    ax[ip].set_ylim(ymin,ymax)
    ax[ip].set_xlim(xticks[ip][0],xticks[ip][-1])
    ax[ip].set_xticks(xticks[ip],xlabels[ip],color=colpal[4])
    ax[ip].tick_params(axis="x", direction="in", length=0)
    ax[ip].tick_params(axis="y", left=False, right=False, direction="in", color=colpal[2], labelcolor=colpal[4], width=1, zorder=0)
    if ip!=0 :
        ax[ip].yaxis.set_ticklabels([])
    for edge in ["bottom", "top", "left", "right"] :
        ax[ip].spines[edge].set_color(colpal[4])
        ax[ip].spines[edge].set_linewidth(1)
        ax[ip].spines[edge].set_zorder(4)
ax[0].tick_params(axis="y", left=True)
ax[-1].tick_params(axis="y", right=True)

f4=open("label.dat","w")
for ip in range(len(xticks)):
    for ik in range(len(xticks[ip])):
        f4.write(str(xticks[ip][ik])+" "+xlabels[ip][ik]+"\n")
f4.close()

fig.savefig(outputname,dpi=1200)

