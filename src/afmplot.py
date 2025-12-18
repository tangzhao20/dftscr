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

ltilt = False
if "tilt" in sys.argv:
    ltilt = True
    sys.argv.remove("tilt")
latom = False
if "atom" in sys.argv:
    latom = True
    sys.argv.remove("atom")
lbohr = False
if "bohr" in sys.argv:
    lbohr = True
    sys.argv.remove("bohr")
lverbose = False
if "verbose" in sys.argv:
    lverbose = True
    sys.argv.remove("verbose")
ltoten = False
if "toten" in sys.argv:
    ltoten = True
    sys.argv.remove("toten")
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
niter = 1000  # number of iterations
alpha = 0.2  # damping factor in the iterative solver
h = 0.2  # step size in the finite difference method, in units of A

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
        elif word[0] == "niter":
            niter = int(word[1])
        elif word[0] == "alpha":
            alpha = float(word[1])
        elif word[0] == "h":
            h = float(word[1])
        elif word[0] in ["fdet", "boundary", "spin"]:
            pass
        else:
            print("Warning: keyword \""+word[0]+"\" is not defined.")

k_spring = k_spring * angstrom**2/electron  # convert spring constant to eV/A^2
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
if icenter < 1 or icenter > nz-2:
    print(f"Error: z-plane index must be in [2, {nz-1}].")
    sys.exit()

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

# ==================== caclulate forces for tilt correction ====================
if ltilt:
    # interpolation from a 2d grid
    toten_2d = []
    for iz in range(3):
        toten_2d.append(scipy.interpolate.RectBivariateSpline(y_grid, x_grid, toten[icenter+iz-1]))  # bicubic spline

    # We use an iterative method to solve D=F(D)/k. For conventional D=F/k, set Niter=1 and alpha=1
    x_init = x_grid[np.newaxis, :]
    y_init = y_grid[:, np.newaxis]
    x_new = np.tile(x_init, (ny, 1))
    y_new = np.tile(y_init, (1, nx))
    for iiter in range(niter):
        fx = (x_init - x_new) * k_spring
        fy = (y_init - y_new) * k_spring
        fx += (toten_2d[1](y_new, x_new-h*0.5, grid=False) - toten_2d[1](y_new, x_new+h*0.5, grid=False)) / h
        fy += (toten_2d[1](y_new-h*0.5, x_new, grid=False) - toten_2d[1](y_new+h*0.5, x_new, grid=False)) / h
        x_incr = fx / k_spring
        y_incr = fy / k_spring
        x_new += x_incr * alpha
        y_new += y_incr * alpha
        if np.max(x_incr**2 + y_incr**2) < 1.e-8:
            if lverbose:
                print(f" tilt correction converge at iiter = {iiter:0d}")
            break
        if iiter == niter-1:
            print(f" tilt correction does not converge in niter = {niter:0d}")

# ==================== calculate kts ====================
if ltilt:
    kts = (toten_2d[0](y_new, x_new, grid=False) - 2*toten_2d[1](y_new, x_new, grid=False)
           + toten_2d[2](y_new, x_new, grid=False)) / z_spacing**2
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

    zmax = apc[:, 2].max()

    atom_x = []
    atom_y = []
    atom_color = []
    edge_color = []
    for ia in range(len(apc)):
        if apc[ia, 2] > zmax-1.0:
            atom_x.append(apc[ia, 0]/funit)
            atom_y.append(apc[ia, 1]/funit)
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
im = ax0.imshow(kts, interpolation='spline36', cmap="Blues_r",
                origin="lower", extent=im_extent, aspect='equal', zorder=1)

if latom:
    ax0.scatter(atom_x, atom_y, c=atom_color, s=12, edgecolors=edge_color, linewidths=0.25, zorder=3)
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
    gs1 = fig1.add_gridspec(1, 1, left=0.14, right=0.7526, top=0.95, bottom=0.15)
    ax2 = gs1.subplots()

    q = ax2.quiver(x_grid/funit, y_grid/funit, (x_new-x_init)/funit, (y_new-y_init)/funit, angles='xy',
                   scale_units='xy', scale=1, color=palette["darkblue"], linewidth=1, zorder=3)
    if lverbose:
        print("iz: "+str(icenter+1)+"    Delta_max: "+f"{np.max((x_new-x_init)**2+(y_new-y_init)**2)**0.5:0.4f} Å")

    if latom:
        ax2.scatter(atom_x, atom_y, c=atom_color, s=24, edgecolors=palette["black"], linewidths=0.25, zorder=2)
    ax2.set_xlim([x_range[0]/funit, x_range[1]/funit])
    ax2.set_ylim([y_range[0]/funit, y_range[1]/funit])
    ax2.set_aspect('equal')

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


# ==================== toten images ====================
if ltoten:
    fig2 = plt.figure(figsize=(5, 3.75))
    gs2 = fig2.add_gridspec(1, 2, wspace=0.02, hspace=0.00, left=0.14, right=0.80,
                            top=0.95, bottom=0.15, width_ratios=[0.6, 0.04])
    [ax3, ax4] = gs2.subplots()

    im_extent = [x_range[0]-x_spacing*0.5, x_range[1]+x_spacing*0.5, y_range[0]-y_spacing*0.5, y_range[1]+y_spacing*0.5]
    for ic in range(len(im_extent)):
        im_extent[ic] = im_extent[ic]/funit
    im = ax3.imshow(toten[icenter], interpolation='nearest', cmap="YlOrBr_r",
                    origin="lower", extent=im_extent, aspect='equal', zorder=1)

    if latom:
        ax3.scatter(atom_x, atom_y, c=atom_color, s=12, edgecolors=edge_color, linewidths=0.25, zorder=3)
    ax3.set_xlim([x_range[0]/funit, x_range[1]/funit])
    ax3.set_ylim([y_range[0]/funit, y_range[1]/funit])
    if lbohr:
        ax3.set_xlabel(r"$\mathit{x}\ (Bohr)$", color=palette["black"])
        ax3.set_ylabel(r"$\mathit{y}\ (Bohr)$", color=palette["black"])
    else:
        ax3.set_xlabel(r"$\mathit{x}\ (Å)$", color=palette["black"])
        ax3.set_ylabel(r"$\mathit{y}\ (Å)$", color=palette["black"])

    cb2 = fig2.colorbar(im, cax=ax4, orientation='vertical')
    cb2.outline.set_linewidth(1)
    cb2.outline.set_color(palette["black"])
    if lbohr:
        ax4.set_ylabel(r"$\mathit{k}_{ts}\ (a.u.)$", color=palette["black"])
    else:
        ax4.set_ylabel("Energy (eV)", color=palette["black"])

    ax3.tick_params(axis="x", bottom=True, right=False, direction="in",
                    color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
    ax3.tick_params(axis="y", left=True, right=False, direction="in",
                    color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
    ax4.tick_params(axis="x", bottom=False, right=False, direction="in",
                    color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
    ax4.tick_params(axis="y", left=False, right=True, direction="in",
                    color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
    for edge in ["bottom", "top", "left", "right"]:
        ax3.spines[edge].set_color(palette["black"])
        ax3.spines[edge].set_linewidth(1)
        ax3.spines[edge].set_zorder(4)

    filename = "toten"
    if latom:
        filename += "_atom"
    if lbohr:
        filename += "_bohr"
    filename += "_"+str(icenter+1)
    filename += ".png"
    fig2.savefig(filename, dpi=1200)
