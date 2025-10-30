#!/usr/bin/env python3

# mae_bands.py med easy (E1) (E2)
# Making the band structure plot with mae decomposition

# Only for VASP outputs

# if E2 exists, the energy range is (E1,E2)
# if only E1 exists, the energy range is (-E1,E1)
# if neither exists, the energy range is (-5 eV, 5 eV)

import sys
import numpy as np
from classes import *

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

x_ticks = kpoints0.xticks_out(rlc)
x_labels = kpoints0.xlabels_out()

x_str = ["x", "y", "z"]

i_med = x_str.index(sys.argv[1])
i_easy = x_str.index(sys.argv[2])

y_range = [-5.0, 5.0]
if len(sys.argv) >= 5:
    y_range[0] = min(float(sys.argv[3]), float(sys.argv[4]))
    y_range[1] = max(float(sys.argv[3]), float(sys.argv[4]))
elif len(sys.argv) == 4:
    y_range[0] = -float(sys.argv[3])
    y_range[1] = float(sys.argv[3])

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

bands = BandsPlot(x_ticks, x_labels, y_range)
bands.plot_bands(eigenval0, kpoints0, rlc)

outputname = "mae_bs.png"

# if eigenval1.Ns==2 :
#    plt.legend()

colors = np.where(mae_proj > 0, bands.palette["red"], bands.palette["green"])
dot_size = 3e3
proj_plot_size = (mae_proj * dot_size)**2
size_thereshold = 0.3
proj_plot_size[proj_plot_size < size_thereshold] = 0

for ispin in range(eigenval0.Ns):
    bands.add_scatter(x, energy[ispin], proj_plot_size[ispin], color=colors[ispin], zorder=3.5-ispin)

bands.fig.savefig(outputname, dpi=1200)
