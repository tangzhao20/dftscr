#!/usr/bin/env python3

# mae_bands.py med easy (E1) (E2)
# Making the band structure plot with mae decomposition

#

# if E2 exists, the energy range is (E1,E2)
# if only E1 exists, the energy range is (-E1,E1)
# if neither exists, the energy range is (-5 eV, 5 eV)

import sys
import numpy as np
from classes import *
from load_data import load_palette
import matplotlib.pyplot as plt
import matplotlib as mpl

poscar0 = Poscar()
eigenval0 = Eigenval()
kpoints0 = KpointsBand()

# Input: EIGENVAL, KPOINTS, POSCAR, (DOSCAR)

poscar0.read_vasp()
rlc = poscar0.rlc()

eigenval0.read_vasp()

kpoints0.read_vasp()

procar0 = Procar()
procar0.read_vasp()

if eigenval0.is_semic == True:
    eigenval0.eigshift(eigenval0.vbm)
    eigenval0.writegap(kpoints0)
else:
    doscar0 = Doscar()
    doscar0.read_vasp()
    eigenval0.eigshift(doscar0.ef)

x = eigenval0.eig_x(kp=kpoints0, rlc=rlc)
energy = eigenval0.eigtrans()

xticks = kpoints0.xticks_out(rlc)
xlabels = kpoints0.xlabels_out()


# width is the physical width of each plots
# Nx is the number of x points in each plots
# ixl and ixr are the left and right index of each panel, to seperate the bands
width = []
Nx = []
ixl = []
ixr = []
for p in x:
    width.append(max(p))
    Nx.append(len(p))
    if ixl == []:
        ixl = [0]
    else:
        ixl.append(ixr[-1])
    ixr.append(ixl[-1]+Nx[-1])

x_str = ["x", "y", "z"]

i_med = x_str.index(sys.argv[1])
i_easy = x_str.index(sys.argv[2])

if len(sys.argv) >= 5:
    ymax = float(sys.argv[4])
    ymin = float(sys.argv[4])
elif len(sys.argv) == 4:
    ymax = float(sys.argv[3])
    ymin = -float(sys.argv[3])
else:
    ymax = 5.0
    ymin = -5.0

# calculate MAE decomposition here
atom_list = poscar0.atom_list()

# SOC factors taken from [New J. Phys. 21, 73054 (2019)]
soc_factor_dict = {"Fe": 0.05965, "Co": 0.07412}

atom_mask = np.isin(atom_list, list(soc_factor_dict))
soc_factors = np.array([soc_factor_dict[a] for a in np.array(atom_list)[atom_mask]])

Ld = np.zeros((3, 5, 5))  # i: x, y, z
# Order of orbs: dxy, dyz, dz2, dxz, x2-y2

Ld[0, 0, 3] = 1.0
Ld[0, 1, 2] = 3.0**0.5
Ld[0, 1, 4] = 1.0
Ld[0, 2, 1] = -3.0**0.5
Ld[0, 3, 0] = -1.0
Ld[0, 4, 1] = -1.0

Ld[1, 0, 1] = 1.0
Ld[1, 1, 0] = -1.0
Ld[1, 2, 3] = -3.0**0.5
Ld[1, 3, 2] = 3.0**0.5
Ld[1, 3, 4] = -1.0
Ld[1, 4, 3] = 1.0

Ld[2, 0, 4] = -2.0
Ld[2, 1, 3] = -1.0
Ld[2, 3, 1] = 1.0
Ld[2, 4, 0] = 2.0

L_med = Ld[i_med, :, :]
L_easy = Ld[i_easy, :, :]

eta = 0.0001

mae_proj = np.zeros((2, eigenval0.Nb, eigenval0.Nk))

for ispin1 in range(2):
    for ispin2 in range(2):
        if ispin1 == ispin2:
            sign = 1
        else:
            sign = -1
        # vectorization over bands [ik, ib1, ib2]
        e_diff = procar0.eig[ispin1, :, :][:, :, np.newaxis] - procar0.eig[ispin2, :, :][:, np.newaxis, :]
        e = e_diff / (e_diff**2 + eta**2)
        f = procar0.occ[ispin1, :, :][:, :, np.newaxis] * (1 - procar0.occ[ispin2, :, :][:, np.newaxis, :])
        # f = procar0.occ[ispin1, ik, :][:, np.newaxis] - procar0.occ[ispin2, ik, :][np.newaxis, :]

        c1 = procar0.complex[ispin1, :, :, atom_mask, 4:9].conj()  # numpy move the masked axis to front
        c2 = procar0.complex[ispin2, :, :, atom_mask, 4:9]
        # a: atoms; i,j: bands; x: directions; m,n: orbitals
        L = np.abs(np.einsum("a, akim, mn, akjn -> kij", soc_factors, c1, L_med, c2))**2 - \
            np.abs(np.einsum("a, akim, mn, akjn -> kij", soc_factors, c1, L_easy, c2))**2
        mae_proj[ispin1, :, :] += sign * np.sum(f * e * L,  axis=2).T * 0.25
        mae_proj[ispin2, :, :] += sign * np.sum(f * e * L,  axis=1).T * 0.25
mae_proj *= 0.5

# make the plot

palette = load_palette()
mpl.rcParams["font.sans-serif"].insert(0, "Noto Sans")
mpl.rcParams.update({'font.size': 14})

spin_label = ["spin up", "spin down"]
line_color = ["darkblue", "orange"]

fig = plt.figure(figsize=(5, 3.75))
gs0 = fig.add_gridspec(1, len(xticks), wspace=0.0, hspace=0.00, left=0.14, right=0.98,
                       top=0.97, bottom=0.07, width_ratios=width[:len(xticks)])
ax = []
for ip in range(len(xticks)):
    ax.append(fig.add_subplot(gs0[ip]))

    ax[ip].grid(axis="x", linewidth=1, color=palette["gray"], zorder=0)
    ax[ip].axhline(linewidth=1, color=palette["gray"], zorder=0)
    for ispin in range(eigenval0.Ns):
        for ib in range(eigenval0.Nb):
            if eigenval0.Ns == 2 and ib == 0:
                ax[ip].plot(x[ip], energy[ispin][ib][ixl[ip]:ixr[ip]], color=palette[line_color[ispin]],
                            label=spin_label[ispin], linewidth=1, zorder=3-ispin)
            else:
                ax[ip].plot(x[ip], energy[ispin][ib][ixl[ip]:ixr[ip]], color=palette[line_color[ispin]],
                            linewidth=1, zorder=3-ispin)

outputname = "mae_bs.png"

# if eigenval1.Ns==2 :
#    plt.legend()

# projection plot
proj_color = ["blue", "beige"]

dot_size = 4e3
proj_plot_size = mae_proj * dot_size
size_thereshold = 0.3
proj_plot_size[np.abs(proj_plot_size) < size_thereshold] = 0

for ispin in range(eigenval0.Ns):
    for ip in range(len(xticks)):
        for ib in range(eigenval0.Nb):
            colors = np.where(mae_proj[ispin, ib, ixl[ip]:ixr[ip]] > 0, palette["red"], palette["green"])
            ax[ip].scatter(x[ip], energy[ispin][ib][ixl[ip]:ixr[ip]],
                           s=np.abs(proj_plot_size[ispin, ib, ixl[ip]:ixr[ip]]),
                           c=colors, zorder=3.5-ispin)


ax[0].set_ylabel("Energy (eV)", labelpad=-2, color=palette["black"])
for ip in range(len(xticks)):
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

fig.savefig(outputname, dpi=1200)
