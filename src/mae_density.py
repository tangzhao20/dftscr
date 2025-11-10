#!/usr/bin/env python3

# mae_bands.py med easy (E1) (E2)
# Plot energy-resolved MAE density

# Only for VASP outputs

# if E2 exists, the energy range is [E1, E2]
# if only E1 exists, the energy range is [-E1, E1]
# if neither exists, the energy range is [-5 eV, 5 eV]

import sys
import numpy as np
from classes import Poscar, Procar, Doscar
from load_data import load_palette

import matplotlib as mpl
import matplotlib.pyplot as plt

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

poscar0 = Poscar()
poscar0.read_vasp(filename="CONTCAR")

procar0 = Procar()
procar0.read_vasp()

doscar0 = Doscar()  # for Fermi level
doscar0.read_vasp()
procar0.eig = procar0.eig - doscar0.ef

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

Ne = 1000
d_energy = np.linspace(y_range[0], y_range[1], Ne)
mae_d = np.zeros((2, Ne))  # spin, energy

sigma = 0.05  # eV
gaussian_coeff = (1 / (sigma * (2*np.pi)**0.5))

for ik in range(procar0.Nk):
    weight = procar0.weight[ik]
    for ispin1 in range(2):
        for ispin2 in range(2):
            if ispin1 == ispin2:
                sign = 1
            else:
                sign = -1

            # vectorization over bands
            e1 = procar0.eig[ispin1, ik, :]
            e2 = procar0.eig[ispin2, ik, :]
            e_diff = e1[:, np.newaxis] - e2[np.newaxis, :]
            e = e_diff / (e_diff**2 + eta**2)
            f = procar0.occ[ispin1, ik, :][:, np.newaxis] * (1 - procar0.occ[ispin2, ik, :][np.newaxis, :])
            # f = procar0.occ[ispin1, ik, :][:, np.newaxis] - procar0.occ[ispin2, ik, :][np.newaxis, :]

            c1 = procar0.complex[ispin1, ik, :, atom_mask, 4:9].conj()  # numpy move the masked axis to front
            c2 = procar0.complex[ispin2, ik, :, atom_mask, 4:9]

            # L matrix [Nb, Nb]
            L = np.abs(np.einsum("a, aim, mn, ajn -> ij", soc_factors, c1, L_med, c2))**2 - \
                np.abs(np.einsum("a, aim, mn, ajn -> ij", soc_factors, c1, L_easy, c2))**2

            d1 = (d_energy[:, None] - e1[None, :]) / sigma  # [Ne, Nb]
            d2 = (d_energy[:, None] - e2[None, :]) / sigma  # [Ne, Nb]
            s1 = gaussian_coeff * np.exp(-0.5 * d1**2)
            s2 = gaussian_coeff * np.exp(-0.5 * d2**2)

            mae_matrix = f * e * L  # [Nb, Nb]

            mae_d[ispin1] += sign * weight * 0.25 * np.einsum("ij, ei -> e", mae_matrix, s1)
            mae_d[ispin2] += sign * weight * 0.25 * np.einsum("ij, ej -> e", mae_matrix, s2)

mae_d = mae_d * 0.5

# print(np.sum(mae_d*(d_energy[1]-d_energy[0])))  # check the integration gives the correct MAE

mae_d_max = np.max(mae_d)
mae_d_min = np.min(mae_d)

palette = load_palette()
mpl.rcParams["font.sans-serif"].insert(0, "Noto Sans")
mpl.rcParams.update({'font.size': 14})

# increase the top space for title
fig = plt.figure(figsize=(1, 3.9))
gs0 = fig.add_gridspec(1, 1, wspace=0.0, hspace=0.00, left=0.03, right=0.97,
                       top=0.9326923076923080, bottom=0.0673076923076923)

ax = fig.add_subplot(gs0[0])

ax.axhline(linewidth=1, color=palette["gray"], zorder=0)
ax.axvline(linewidth=1, color=palette["gray"], zorder=0)

ax.plot(mae_d[0], d_energy, color=palette["darkblue"], linewidth=1, zorder=3.5)
ax.plot(mae_d[1], d_energy, color=palette["orange"], linewidth=1, zorder=3)

x_range = [0.0, 0.0]
x_range[0] = mae_d_min - (mae_d_max - mae_d_min) * 0.05
x_range[1] = mae_d_max + (mae_d_max - mae_d_min) * 0.05
ax.set_xlim(x_range)

ax.set_ylim(y_range)
ax.tick_params(axis="x", bottom=False, top=False, direction="in", length=0)
ax.tick_params(axis="y", left=True, right=True, direction="in", color=palette["gray"],
               labelcolor=palette["black"], width=1, zorder=0, pad=4)
ax.set_yticklabels([])

ax.set_title("MAE den.", fontsize=14, pad=4, color=palette["black"])
ax.set_xticks([mae_d_min, 0, mae_d_max], ["–", "0", "+"], color=palette["black"])

for edge in ["bottom", "top", "left", "right"]:
    ax.spines[edge].set_color(palette["black"])
    ax.spines[edge].set_linewidth(1)
    ax.spines[edge].set_zorder(4)

fig.savefig("mae_d.png", dpi=1200)
