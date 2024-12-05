#!/usr/bin/env python3

# Plot the contour of kts vs y and z.

# This script only works for a path on the y axis now.
# In the afm.in, set x_range to two same numbers, for example:
# x_range 0.0 0.0

import sys
import os
from load_data import load_constant, load_palette, load_atom_color
from classes import Poscar
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
import numpy as np

bohr = load_constant("bohr")
rydberg = load_constant("rydberg")
Ha = rydberg*2.0
electron = load_constant("electron")
angstrom = load_constant("angstrom")

latom = False
if "atom" in sys.argv:
    latom = True
    sys.argv.remove("atom")
lmax = False
if "max" in sys.argv:
    lmax = True
    sys.argv.remove("max")
lbohr = False
if "bohr" in sys.argv:
    lbohr = True
    sys.argv.remove("bohr")

# ==================== read the input file ====================
x_spacing = 0.6
y_spacing = 0.6
z_spacing = 0.3
z_range = [5.7, 6.3]
parallel = 1
k_spring = 0.8  # k in N/m
k_spring = k_spring * angstrom**2/electron  # convert spring constant to eV/A^2

with open("afm.in", "r") as f0:
    for l in f0:
        word = l.split()
        if not word or word[0][0] in ("#", "!"):
            continue
        if word[0] == "x_range":
            x_range = [float(word[1]), float(word[2])]
        elif word[0] == "y_range":
            y_range = [float(word[1]), float(word[2])]
        elif word[0] == "z_range":
            z_range = [float(word[1]), float(word[2])]
        elif word[0] == "x_spacing":
            x_spacing = float(word[1])
        elif word[0] == "y_spacing":
            y_spacing = float(word[1])
        elif word[0] == "z_spacing":
            z_spacing = float(word[1])
        elif word[0] == "parallel":
            parallel = int(word[1])
        elif word[0] == "k_spring":
            k_spring = float(word[1])
        elif word[0] in ["fdet", "boundary", "spin"]:
            pass
        else:
            print("Warning: keyword \""+word[0]+"\" is not defined.")

x_spacing = x_spacing * bohr
y_spacing = y_spacing * bohr
z_spacing = z_spacing * bohr
x_range = [x_range[0] * bohr, x_range[1] * bohr]
y_range = [y_range[0] * bohr, y_range[1] * bohr]
z_range = [z_range[0] * bohr, z_range[1] * bohr]

# ==================== prepare the x and y coordinates ====================
x = x_range[0]
x_grid = np.arange(x_range[0], x_range[1]+1e-6, x_spacing)
nx = len(x_grid)

y = y_range[0]
y_grid = np.arange(y_range[0], y_range[1]+1e-6, y_spacing)
ny = len(y_grid)

z = z_range[0]
z_grid = np.arange(z_range[0], z_range[1]+1e-6, z_spacing)
nz = len(z_grid)

# ==================== calculate or read toten ====================
toten = np.zeros((nz, ny, nx))  # in eV
files = os.listdir()
# if the toten.dat exist, read it, if not, read from calculations outputs
if "toten.dat" in files:
    il = 0
    with open("toten.dat", "r") as f2:
        for l in f2:
            word = l.split()
            if not word or word[0][0] in ("#", "!"):
                continue
            ix, iy, iz = map(int, word[:3])
            toten[iz, iy, ix] = float(word[3])
            il += 1
    if il != nx * ny * nz:
        print(f"Warning: Expected {nx * ny * nz} data points, but found {il} in toten.dat")

else:
    # Set up AFM scan path
    step_list = []
    with open("steps.dat", "r") as f3:
        for l in f3:
            word = l.split()
            if not word or word[0][0] in ("#", "!"):
                continue
            step_list.append(int(word[0]))

    scan_path = np.zeros((ny*nx, 2), int)  # index [iy, ix] of each point
    scan_path[:, 0] = np.repeat(range(ny), nx)
    scan_path[:, 1] = np.tile(np.concatenate([np.arange(nx), np.arange(nx-1, -1, -1)]), ny//2+1)[:ny*nx]

    for iz in range(nz):
        istep = -1
        for ip in range(parallel):
            filename = "seq_"+str(iz+1)+"_"+str(ip+1)+"/parsec.out"
            f1 = open(filename, "r")
            line = f1.readlines()
            f1.close()
            with open("seq_"+str(iz+1)+"_"+str(ip+1)+"/parsec.out", "r") as f1:
                for l in f1:
                    word = l.split()
                    if not word or word[0][0] in ("#", "!"):
                        continue
                    if len(word) >= 2 and word[0] == "Starting" and word[1] == "SCF...":
                        istep += 1
                    if len(word) >= 5 and word[0] == "Total" and word[1] == "Energy" and word[2] == "=":
                        toten[iz, scan_path[istep][0], scan_path[istep][1]] = float(word[3]) * rydberg

    # write toten.dat
    with open("toten.dat", "w") as f2:
        f2.write("#   ix    iy    iz      toten(eV)\n")
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    f2.write(f"{ix:6d}{iy:6d}{iz:6d}{toten[iz][iy][ix]:24.12f}\n")

# ==================== calculate kts ====================
kts = (toten[:nz-2, :, 0] - 2*toten[1:nz-1, :, 0] + toten[2:nz, :, 0]) / z_spacing**2
if lbohr:
    # convert k_ts from eV/A^2 to Ha/a0^2
    kts = kts * bohr**2/Ha

# ==================== calculate maximums ====================
if lmax:
    interp_func = scipy.interpolate.RectBivariateSpline(z_grid[1:nz-1], y_grid, kts, kx=3, ky=3)
    z_max_list = np.linspace(z_range[0], z_range[1], 500)
    y_max_list = []
    for z in z_max_list:
        def y_grid_1d(y): return -interp_func(z, y)[0]  # Negative for maximizing
        z_max_result = scipy.optimize.minimize_scalar(y_grid_1d, bounds=(y_range[0], y_range[1]), method='bounded')
        y_max_list.append(z_max_result.x[0])

# ==================== construct the atomic structure ====================
palette = load_palette()
if lbohr:
    funit = bohr
else:
    funit = 1.0
if latom:
    color_dict = load_atom_color()
    poscar1 = Poscar()
    poscar1.read_parsec(filename="sample.parsec_st.dat")
    # poscar1 is modified here. If a feature needs to keep poscar1, use copy.
    if poscar1.Ndim == 2:
        poscar1.supercell([2, 2, 1])
        poscar1.move([-0.5, -0.5, 0], lbox=False)
    apc = poscar1.cartesian()
    atom = poscar1.atom_list()

    zmax = -1e6
    for ia in range(len(apc)):
        zmax = max(zmax, apc[ia][2])

    atom_y = []
    atom_color = []
    for ia in range(len(apc)):
        if apc[ia][2] > zmax - 1.0 and abs(apc[ia][0]) < 1.e-4:
            atom_y.append(apc[ia][1]/funit)
            atom_color.append(palette[color_dict[atom[ia]]])

# ==================== creating the plot ====================
mpl.rcParams["font.sans-serif"].insert(0, "Noto Sans")
mpl.rcParams.update({'font.size': 14})
mpl.rcParams.update({'mathtext.default': 'regular'})

fig0 = plt.figure(figsize=(5, 3.75))
gs0 = fig0.add_gridspec(1, 2, wspace=0.02, hspace=0.00, left=0.14, right=0.80,
                        top=0.95, bottom=0.15, width_ratios=[0.6, 0.04])
[ax0, ax1] = gs0.subplots()

im_extent = [y_range[0] - y_spacing*0.5, y_range[1] + y_spacing *
             0.5, z_range[0] - z_spacing*0.5, z_range[1] + z_spacing*0.5]
for ic in range(len(im_extent)):
    im_extent[ic] = im_extent[ic]/funit
im = ax0.imshow(kts, interpolation='bicubic', cmap="rainbow", origin="lower", extent=im_extent, zorder=1)
ax0.set_aspect((y_range[1]-y_range[0])/(z_range[1]-z_range[0]))

if latom:
    for ia in range(len(atom_y)):
        ax0.axvline(x=atom_y[ia], color=atom_color[ia], linestyle="dashdot", linewidth=1)
if lmax:
    ax0.plot(y_max_list, z_max_list, color=palette["red"], linestyle="dotted", linewidth=1)

ax0.set_xlim([y_range[0]/funit, y_range[1]/funit])
ax0.set_ylim([z_range[0]/funit, z_range[1]/funit])
if lbohr:
    ax0.set_xlabel(r"$\mathit{y}\ (Bohr)$", color=palette["black"])
    ax0.set_ylabel(r"$\mathit{z}\ (Bohr)$", color=palette["black"])
else:
    ax0.set_xlabel(r"$\mathit{y}\ (Å)$", color=palette["black"])
    ax0.set_ylabel(r"$\mathit{z}\ (Å)$", color=palette["black"])

cb = fig0.colorbar(im, cax=ax1, orientation='vertical')
cb.outline.set_linewidth(1)
cb.outline.set_color(palette["black"])
if lbohr:
    ax1.set_ylabel(r"$\mathit{k}_{ts}\ (a.u.)$", color=palette["black"])
else:
    ax1.set_ylabel(r"$\mathit{k}_{ts}\ (eV/Å^2)$", color=palette["black"])

ax0.tick_params(axis="x", bottom=True, right=False, direction="in",
                color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
ax0.tick_params(axis="y", left=True, right=False, direction="in",
                color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
ax1.tick_params(axis="x", bottom=False, right=False, direction="in",
                color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
ax1.tick_params(axis="y", left=False, right=True, direction="in",
                color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
for edge in ["bottom", "top", "left", "right"]:
    ax0.spines[edge].set_color(palette["black"])
    ax0.spines[edge].set_linewidth(1)
    ax0.spines[edge].set_zorder(4)

filename = "afm_line"
if latom:
    filename += "_atom"
if lmax:
    filename += "_max"
if lbohr:
    filename += "_bohr"
filename += ".png"
fig0.savefig(filename, dpi=1200)
