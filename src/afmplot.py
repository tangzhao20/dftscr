#!/usr/bin/env python3

# This script reads the afm simulation results and makes the plot

# The matrix in this script are mostly M[ny][nx] or M[nz][ny][nx] to match the size of image

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
lbohr = False
if "bohr" in sys.argv:
    lbohr = True
    sys.argv.remove("bohr")
ltilt = False
if "tilt" in sys.argv:
    ltilt = True
    sys.argv.remove("tilt")
icenter = 1
for word in sys.argv:
    if word.isnumeric():
        icenter = int(word)-1
        sys.argv.remove(word)

# ==================== read the input file ====================
x_spacing = 0.6
y_spacing = 0.6
z_spacing = 0.3
z_range = [5.7, 6.3]
parallel = 1
k_spring = 0.8  # k in N/m
k_spring = k_spring * angstrom**2/electron  # convert spring constant to eV/A^2

f0 = open("afm.in", "r")
line = f0.readlines()
f0.close()
for l in line:
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
    f3 = open("steps.dat", "r")
    line = f3.readlines()
    f3.close()
    step_list = []
    for l in line:
        word = l.split()
        if not word or word[0][0] in ("#", "!"):
            continue
        step_list.append(int(word[0]))
    lxincrease = True
    k = 0
    ip = 0
    scan_path = [] # index (ix,iy) of each point
    for iy in range(ny):
        if lxincrease:
            for ix in range(nx):
                if k == 0:
                    scan_path.append([])
                scan_path[-1].append([ix, iy])
                k += 1
                if k == step_list[ip]:
                    ip += 1
                    k = 0
        else:
            for ix in range(nx-1, -1, -1):
                if k == 0:
                    scan_path.append([])
                scan_path[-1].append([ix, iy])
                k += 1
                if k == step_list[ip]:
                    ip += 1
                    k = 0
        lxincrease = not lxincrease

    for iz in range(nz):
        for ip in range(parallel):
            istep = -1
            filename = "seq_"+str(iz+1)+"_"+str(ip+1)+"/parsec.out"
            f1 = open(filename, "r")
            line = f1.readlines()
            f1.close()

            for l in line:
                word = l.split()
                if not word or word[0][0] in ("#", "!"):
                    continue
                if len(word) >= 2 and word[0] == "Starting" and word[1] == "SCF...":
                    istep += 1
                if len(word) >= 5 and word[0] == "Total" and word[1] == "Energy" and word[2] == "=":
                    toten[iz, scan_path[ip][istep][1], scan_path[ip][istep][0]] = float(word[3]) * rydberg

    # write the toten file
    f2 = open("toten.dat", "w")
    f2.write("#   ix    iy    iz        toten(eV)\n")
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                f2.write(f"{ix:6d}{iy:6d}{iz:6d}{toten[iz][iy][ix]:24.12f}\n")
    f2.close()

# ==================== caclulate forces for tilt correction ====================
if ltilt:
    fx = np.zeros((ny, nx))
    fy = np.zeros((ny, nx))
    for iy in range(ny):
        for ix in range(nx):
            if ix != 0 and ix != nx-1:
                fx[iy][ix] = (toten[icenter][iy][ix-1]-toten[icenter][iy][ix+1])*0.5/x_spacing
            elif ix == 0:
                fx[iy][ix] = (toten[1][iy][ix]-toten[1][iy][ix+1])/x_spacing
            elif ix == nx-1:
                fx[iy][ix] = (toten[1][iy][ix-1]-toten[1][iy][ix])/x_spacing
            if iy != 0 and iy != ny-1:
                fy[iy][ix] = (toten[icenter][iy-1][ix]-toten[icenter][iy+1][ix])*0.5/y_spacing
            elif iy == 0:
                fy[iy][ix] = (toten[1][iy][ix]-toten[1][iy+1][ix])/y_spacing
            elif iy == ny-1:
                fy[iy][ix] = (toten[1][iy-1][ix]-toten[1][iy][ix])/y_spacing

    x_new = np.zeros((ny, nx))
    y_new = np.zeros((ny, nx))
    for iy in range(ny):
        for ix in range(nx):
            x_new[iy, ix] = x_grid[ix] + fx[iy][ix] / k_spring
            y_new[iy, ix] = y_grid[iy] + fy[iy][ix] / k_spring

    # create 2d map from toten
    toten_2d = []
    for iz in range(nz):
        toten_2d0 = scipy.interpolate.RectBivariateSpline(y_grid, x_grid, toten[iz])
        toten_2d.append(toten_2d0)

# ==================== calculate kts ====================
kts = np.zeros((ny, nx))
if ltilt:
    for iy in range(ny):
        for ix in range(nx):
            kts[iy][ix] = (toten_2d[icenter-1](y_new[iy, ix], x_new[iy, ix])[0, 0]
                           - 2*toten_2d[icenter](y_new[iy, ix], x_new[iy, ix])[0, 0]
                           + toten_2d[icenter+1](y_new[iy, ix], x_new[iy, ix])[0, 0]) / z_spacing**2
else:
    kts = (toten[icenter-1, :, :] - 2*toten[icenter, :, :] + toten[icenter+1, :, :]) / z_spacing**2

if lbohr:
    # convert k_ts from eV/A^2 to Ha/a0^2
    kts = kts * bohr**2/Ha

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

    atom_x = []
    atom_y = []
    atom_color = []
    edge_color = []
    for ia in range(len(apc)):
        if apc[ia][2] > zmax-1.0:
            atom_x.append(apc[ia][0]/funit)
            atom_y.append(apc[ia][1]/funit)
            atom_color.append(palette[color_dict[atom[ia]]])
            if color_dict[atom[ia]].startswith("dark") or color_dict[atom[ia]] == "black":
                edge_color.append(palette["white"])
            else:
                edge_color.append(palette["black"])

# ==================== creating the plot ====================
mpl.rcParams["font.sans-serif"].insert(0, "Noto Sans")
mpl.rcParams.update({'font.size': 14})
mpl.rcParams.update({'mathtext.default': 'regular'})

fig0 = plt.figure(figsize=(5, 3.75))
gs0 = fig0.add_gridspec(1, 2, wspace=0.02, hspace=0.00, left=0.14, right=0.80,
                        top=0.95, bottom=0.15, width_ratios=[0.6, 0.04])
[ax0, ax1] = gs0.subplots()

im_extent = [x_range[0]-x_spacing*0.5, x_range[1]+x_spacing*0.5, y_range[0]-y_spacing*0.5, y_range[1]+y_spacing*0.5]
for ic in range(len(im_extent)):
    im_extent[ic] = im_extent[ic]/funit
im = ax0.imshow(kts, interpolation='bicubic', cmap="YlOrBr_r",
                origin="lower", extent=im_extent, aspect='equal', zorder=1)

if latom:
    ax0.scatter(atom_x, atom_y, c=atom_color, s=12, edgecolors=edge_color, linewidths=1, zorder=3)
ax0.set_xlim([x_range[0]/funit, x_range[1]/funit])
ax0.set_ylim([y_range[0]/funit, y_range[1]/funit])
if lbohr:
    ax0.set_xlabel(r"$\mathit{x}\ (Bohr)$", color=palette["black"])
    ax0.set_ylabel(r"$\mathit{y}\ (Bohr)$", color=palette["black"])
else:
    ax0.set_xlabel(r"$\mathit{x}\ (Å)$", color=palette["black"])
    ax0.set_ylabel(r"$\mathit{y}\ (Å)$", color=palette["black"])

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

filename = "afm"
if ltilt:
    filename += "_tilt"
if latom:
    filename += "_atom"
if lbohr:
    filename += "_bohr"
filename += "_"+str(icenter+1)
filename += ".png"
fig0.savefig(filename, dpi=1200)

# ==================== vector map for tilt corrections ====================
if ltilt:
    fig1 = plt.figure(figsize=(5, 3.75))
    gs1 = fig1.add_gridspec(1, 1, left=0.14, right=0.74, top=0.95, bottom=0.15)
    ax2 = gs1.subplots()

    q = ax2.quiver(x_grid/funit, y_grid/funit, fx/k_spring/funit, fy/k_spring/funit,
                   color=palette["darkblue"], linewidth=1, zorder=2)

    if latom:
        ax2.scatter(atom_x, atom_y, c=atom_color, s=24, edgecolors="none", linewidths=1, zorder=3)
    ax2.set_xlim([x_range[0]/funit, x_range[1]/funit])
    ax2.set_ylim([y_range[0]/funit, y_range[1]/funit])
    if lbohr:
        ax2.set_xlabel(r"$\mathit{x}\ (Bohr)$", color=palette["black"])
        ax2.set_ylabel(r"$\mathit{y}\ (Bohr)$", color=palette["black"])
    else:
        ax2.set_xlabel(r"$\mathit{x}\ (Å)$", color=palette["black"])
        ax2.set_ylabel(r"$\mathit{y}\ (Å)$", color=palette["black"])

    ax2.tick_params(axis="x", bottom=True, right=False, direction="in",
                    color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
    ax2.tick_params(axis="y", left=True, right=False, direction="in",
                    color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
    for edge in ["bottom", "top", "left", "right"]:
        ax2.spines[edge].set_color(palette["black"])
        ax2.spines[edge].set_linewidth(1)
        ax2.spines[edge].set_zorder(4)

    filename = "tilt"
    if latom:
        filename += "_atom"
    if lbohr:
        filename += "_bohr"
    filename += "_"+str(icenter+1)
    filename += ".png"
    fig1.savefig(filename, dpi=1200)
