#!/usr/bin/env python3

# bands.py package (E1) (E2)
# Making the band structure plot

# if E2 exists, the energy range is (E1,E2)
# if only E1 exists, the energy range is (-E1,E1)
# if neither exists, the energy range is (-5 eV, 5 eV)

import os
import sys
from classes import *
from load_data import load_package_name, load_palette
import matplotlib.pyplot as plt
import matplotlib as mpl

fsecond = False

if len(sys.argv) <= 1:
    print("python bands.py package (Emin) (Emax)")
    sys.exit()
package = sys.argv[1]

package_name = load_package_name()

poscar1 = Poscar()
eigenval1 = Eigenval()
kpoints1 = KpointsBand()

if package in package_name["vasp"]+package_name["vaspproj"]:
    # Input: EIGENVAL, KPOINTS, POSCAR, (DOSCAR)

    poscar1.read_vasp()
    rlc = poscar1.rlc()

    eigenval1.read_vasp()

    kpoints1.read_vasp()

    if eigenval1.is_semic == True:
        eigenval1.eigshift(eigenval1.vbm)
        eigenval1.writegap(kpoints1)
    else:
        doscar1 = Doscar()
        doscar1.read_vasp()
        eigenval1.eigshift(doscar1.ef)

    x = eigenval1.bandkpout(kp=kpoints1, rlc=rlc)
    energy = eigenval1.eigtrans()

    xticks = kpoints1.xticks_out(rlc)
    xlabels = kpoints1.xlabels_out()

elif package in package_name["qe"]+package_name["qeproj"]:
    # Input: *.xml, kpath.in
    # No need to run bands.x

    poscar1.read_xml()
    rlc = poscar1.rlc()

    eigenval1.read_qexml()

    kpoints1.read_kpathin()

    if eigenval1.is_semic == True:
        eigenval1.eigshift(eigenval1.vbm)
        eigenval1.writegap(kpoints1)
    else:
        doscar1 = Doscar()
        doscar1.read_xml()
        eigenval1.eigshift(doscar1.ef)

    x = eigenval1.bandkpout(kp=kpoints1, rlc=rlc)
    energy = eigenval1.eigtrans()

    xticks = kpoints1.xticks_out(rlc)
    xlabels = kpoints1.xlabels_out()

elif package in package_name["wannier90"]:
    # Input : nscf.in, ../bands/*.xml, *_band.kpt, *_band.dat, kpath.in
    # Only semiconductors are supported

    poscar1.read_qe("nscf.in")
    rlc = poscar1.rlc()

    pad = 0
    for w in sys.argv:
        if w.startswith("pad="):
            pad = int(w.split("=")[1])
            sys.argv.remove(w)
            break

    eigenval1.read_wan(Nb_pad=pad)

    kpoints1.read_kpathin()

    # find a ../bands/*.xml file
    files = os.listdir("../bands")
    for f in files:
        if f.endswith('.xml'):
            filename = f
            fsecond = True
            break

    if fsecond:
        eigenval2 = Eigenval()
        eigenval2.read_qexml("../bands/"+filename)

        eigenval2.gap()
        eigenval2.eigshift(eigenval2.vbm)
        if eigenval2.is_semic == True:
            eigenval1.is_semic = True

    eigenval1.occ = []
    for ik in range(eigenval1.Nk):
        occ0 = []
        for ib in range(eigenval2.Nvb[0]-pad):
            occ0.append([1.0])
        if fsecond:
            for ib in range(eigenval2.Nvb[0]-pad, eigenval1.Nb):
                occ0.append([0.0])
            eigenval1.occ.append(occ0)

    if eigenval1.is_semic == True:
        eigenval1.gap()
        eigenval1.writegap(kpoints1)
        eigenval1.eigshift(eigenval1.vbm)
    else:
        print("Metal band structure are not shifted")

    x = eigenval1.bandkpout(kp=kpoints1, rlc=rlc)
    energy = eigenval1.eigtrans()

    xticks = kpoints1.xticks_out(rlc)
    xlabels = kpoints1.xlabels_out()

    if fsecond:
        x2 = eigenval2.bandkpout(kp=kpoints1, rlc=rlc)
        energy2 = eigenval2.eigtrans()


elif package in package_name["parsec"]:
    # Input: bands.dat, parsec.in, kpath.in

    poscar1.read_parsec()
    rlc = poscar1.rlc()

    eigenval1.read_parsec()
    eigenval1.kc2kd(poscar1.lc)

    kpoints1.read_kpathin()

    if eigenval1.is_semic == True:
        eigenval1.gap()
        eigenval1.writegap(kpoints1)
        eigenval1.eigshift(eigenval1.vbm)
    else:
        print("Metal band structure are not shifted")

    x = eigenval1.bandkpout(kp=kpoints1, rlc=rlc)
    energy = eigenval1.eigtrans()

    xticks = kpoints1.xticks_out(rlc)
    xlabels = kpoints1.xlabels_out()

else:
    print("Package \""+package+"\" is not supported yet.")
    print("python bands.py package (Emin) (Emax)")
    sys.exit()

# width is the physical width of each plots
# Nx is the number of x points in each plots
# ixl and ixr are the left and right index of each panel, to seperate the bands
Np = len(xticks)
width = [0.0] * Np
Nx = [0] * Np
ixl = [0] * Np
ixr = [0] * Np
for ip in range(Np):
    width[ip] = xticks[ip][-1]
    Nx[ip] = len(x[ip])
    ixl[ip] = ixr[ip-1]  # ixl[0] = 0
    ixr[ip] = ixl[ip] + Nx[ip]
if fsecond:
    Nx2 = [0] * Np
    ix2l = [0] * Np
    ix2r = [0] * Np
    for p in x2:
        Nx2[ip] = len(x2[ip])
        ix2l[ip] = ix2r[ip-1]  # ix2l[0] = 0
        ix2r[ip] = ix2l[ip] + Nx2[ip]

is_proj = False
if package in package_name["vaspproj"]+package_name["qeproj"]:
    is_proj = True

    procar1 = Procar()
    if package in package_name["vaspproj"]:
        # Input: PROCAR
        procar1.read_vasp()
    elif package in package_name["qeproj"]:
        # Input: projwfc.out
        procar1.read_qe()

    atom_list = sys.argv[2]
    del sys.argv[2]
    atom_flag = poscar1.read_atom_list(atom_list)
    orb_list = sys.argv[2]
    del sys.argv[2]
    orb_flag = procar1.read_orb_list(orb_list)

# The following plotting section should be generalized, not code-specific

if len(sys.argv) >= 4:
    ymax = float(sys.argv[3])
    ymin = float(sys.argv[2])
elif len(sys.argv) == 3:
    ymax = float(sys.argv[2])
    ymin = -float(sys.argv[2])
else:
    ymax = 5.0
    ymin = -5.0

palette = load_palette()
mpl.rcParams["font.sans-serif"].insert(0, "Noto Sans")
mpl.rcParams.update({'font.size': 14})

# band structure plot

spin_label = ["spin up", "spin down"]
line_color = ["darkblue", "orange"]

fig = plt.figure(figsize=(5, 3.75))
gs0 = fig.add_gridspec(1, Np, wspace=0.0, hspace=0.00, left=0.14, right=0.98,
                       top=0.97, bottom=0.07, width_ratios=width[:Np])
ax = []
for ip in range(Np):
    ax.append(fig.add_subplot(gs0[ip]))

    ax[ip].grid(axis="x", linewidth=1, color=palette["gray"], zorder=0)
    ax[ip].axhline(linewidth=1, color=palette["gray"], zorder=0)
    for ispin in range(eigenval1.Ns):
        for ib in range(eigenval1.Nb):
            if eigenval1.Ns == 2 and ib == 0:
                ax[ip].plot(x[ip], energy[ispin][ib][ixl[ip]:ixr[ip]], color=palette[line_color[ispin]],
                            label=spin_label[ispin], linewidth=1, zorder=3-ispin)
            else:
                ax[ip].plot(x[ip], energy[ispin][ib][ixl[ip]:ixr[ip]], color=palette[line_color[ispin]],
                            linewidth=1, zorder=3-ispin)

outputname = "bs.png"

# if eigenval1.Ns==2 :
#    plt.legend()

f3 = open("eigenval.dat", "w")
if eigenval1.Ns == 2:
    for ib in range(len(energy[0])):
        for ip in range(Np):
            ik0 = 0
            for ik in range(Nx[ip]):
                f3.write(str(x[ip][ik])+" "+str(energy[0][ib][ik0])+" "+str(energy[1][ib][ik0])+"\n")
                ik0 += 1
        f3.write("\n")
else:
    for ib in range(len(energy[0])):
        for ip in range(Np):
            ik0 = 0
            for ik in range(Nx[ip]):
                f3.write(str(x[ip][ik])+" "+str(energy[0][ib][ik0])+" "+"\n")
                ik0 += 1
        f3.write("\n")
f3.close()

# projection plot

if is_proj:
    if eigenval1.Ns == 1:
        proj_color = ["orange"]
    else:
        proj_color = ["blue", "beige"]
    dot_size = 50.0
    proj_plot_size = procar1.plot(atom_flag, orb_flag) * dot_size
    for ispin in range(eigenval1.Ns):
        for ip in range(Np):
            for ib in range(eigenval1.Nb):
                ax[ip].scatter(x[ip], energy[ispin][ib][ixl[ip]:ixr[ip]], s=proj_plot_size[ispin, ib, ixl[ip]:ixr[ip]],
                               c=palette[proj_color[ispin]], zorder=2.5-ispin)
    outputname = "proj_"+atom_list+"_"+orb_list+".png"

    f2 = open("proj_"+atom_list+"_"+orb_list+".dat", "w")
    f2.write("#ispin x_kpath energy(eV) size\n")
    for ispin in range(eigenval1.Ns):
        for ib in range(eigenval1.Nb):
            for ip in range(Np):
                ik0 = 0
                for ik in range(Nx[ip]):
                    f2.write(f"{ispin}  {x[ip][ik]}  {energy[ispin][ib][ik0]}  {proj_plot_size[ispin, ib, ik0]}\n")
                    ik0 += 1
            f2.write("\n")
    f2.close()

# second band structure plot (for wannier)

if fsecond:
    for ip in range(Np):
        for ib in range(eigenval2.Nb):
            ax[ip].plot(x2[ip], energy2[0][ib][ix2l[ip]:ix2r[ip]], color=palette["orange"], linewidth=1, zorder=2)
    f3 = open("eigenval2.dat", "w")
    for ib in range(len(energy2[0])):
        for ip in range(Np):
            ik0 = 0
            for ik in range(Nx2[ip]):
                f3.write(str(x2[ip][ik])+" "+str(energy2[0][ib][ik0])+" "+"\n")
                ik0 += 1
        f3.write("\n")
    f3.close()

ax[0].set_ylabel("Energy (eV)", labelpad=-2, color=palette["black"])
for ip in range(Np):
    ax[ip].set_ylim(ymin, ymax)
    ax[ip].set_xlim(xticks[ip][0], xticks[ip][-1])
    ax[ip].set_xticks(xticks[ip], xlabels[ip], color=palette["black"])
    ax[ip].tick_params(axis="x", direction="in", length=0)
    ax[ip].tick_params(axis="y", left=False, right=False, direction="in",
                       color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
    if ip != 0:
        ax[ip].yaxis.set_ticklabels([])
    for edge in ["bottom", "top", "left", "right"]:
        ax[ip].spines[edge].set_color(palette["black"])
        ax[ip].spines[edge].set_linewidth(1)
        ax[ip].spines[edge].set_zorder(4)
ax[0].tick_params(axis="y", left=True)
ax[-1].tick_params(axis="y", right=True)

f4 = open("label.dat", "w")
for ip in range(Np):
    for ik in range(len(xticks[ip])):
        f4.write(str(xticks[ip][ik])+" "+xlabels[ip][ik]+"\n")
f4.close()

fig.savefig(outputname, dpi=1200)
