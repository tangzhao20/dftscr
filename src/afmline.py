#!/usr/bin/env python3

"""
Plot the contour of kts vs y and z.

The current version only works for a path on the y axis.
In the afm.in, set x_range to two same numbers (e.g., x_range 0.0 0.0).

inputs: afm.in, [toten.dat], [steps.dat], [sample.parsec_st.dat]

examples:
    python afmline.py --tilt --atom --max
"""

import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy

from load_data import load_constant, load_palette, load_atom_color
from classes import Poscar


bohr = load_constant("bohr")
rydberg = load_constant("rydberg")
Ha = rydberg*2.0
electron = load_constant("electron")
angstrom = load_constant("angstrom")

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument("--tilt", "-t", action="store_true", help="use tilting correction")
parser.add_argument("--atom", "-a", action="store_true", help="display atom positions")
parser.add_argument("--max", "-m", action="store_true", help="display maximum")
parser.add_argument("--bohr", "-b", action="store_true", help="use atomic units")
parser.add_argument("--force", "-f", action="store_true", help="plot force instead of kts")
parser.add_argument("--verbose", action="store_true", help="print verbose output")

args = parser.parse_args()

use_tilt = args.tilt
show_atom = args.atom
show_max = args.max
use_bohr = args.bohr
plot_force = args.force
verbose = args.verbose

# ==================== read the input file ====================
z_range = [5.7, 6.3]
x_spacing = 0.6
y_spacing = 0.6
z_spacing = 0.3
parallel = 1
k_spring = 0.8  # k in N/m
niter = 1000  # number of iterations
alpha = 0.2  # damping factor in the iterative solver
h = 0.2  # step size in the finite difference method, in units of A

with open("afm.in", "r") as f0:
    for line in f0:
        word = line.split("#")[0].split("!")[0].replace(":", " ").replace("=", " ").split()
        if not word:
            continue

        match word:
            case ["x_range", v1, v2]: x_range = [float(v1), float(v2)]
            case ["y_range", v1, v2]: y_range = [float(v1), float(v2)]
            case ["z_range", v1, v2]: z_range = [float(v1), float(v2)]
            case ["x_spacing", val]: x_spacing = float(val)
            case ["y_spacing", val]: y_spacing = float(val)
            case ["z_spacing", val]: z_spacing = float(val)

            case ["parallel", val]: parallel = int(val)

            case ["k_spring", val]: k_spring = float(val)
            case ["niter", val]: niter = int(val)
            case ["alpha", val]: alpha = float(val)
            case ["h", val]: h = float(val)

            case [("contrast_range" | "fdet" | "boundary" | "spin"), *_]:
                pass

            case _:
                print(f"Warning: check your input line: {' '.join(word)}")

x_range = [x_range[0] * bohr, x_range[1] * bohr]
y_range = [y_range[0] * bohr, y_range[1] * bohr]
z_range = [z_range[0] * bohr, z_range[1] * bohr]
x_spacing = x_spacing * bohr
y_spacing = y_spacing * bohr
z_spacing = z_spacing * bohr
k_spring = k_spring * angstrom**2/electron  # convert spring constant to eV/A^2

# ==================== prepare the x and y coordinates ====================
x = x_range[0]
x_grid = np.arange(x_range[0], x_range[1]+1e-6, x_spacing)
nx = len(x_grid)
if nx > 1:
    print("Warning: nx > 1")

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
        for line in f2:
            word = line.split("#")[0].split("!")[0].split()
            if not word:
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
        for line in f3:
            word = line.split("#")[0].split("!")[0].split()
            if not word:
                continue
            step_list.append(int(word[0]))

    scan_path = np.zeros((ny*nx, 2), int)  # index [iy, ix] of each point
    scan_path[:, 0] = np.repeat(range(ny), nx)
    scan_path[:, 1] = np.tile(np.concatenate([np.arange(nx), np.arange(nx-1, -1, -1)]), ny//2+1)[:ny*nx]

    for iz in range(nz):
        istep = -1
        for ip in range(parallel):
            with open(f"seq_{iz+1}_{ip+1}/parsec.out", "r") as f1:
                for line in f1:
                    if "Starting SCF..." in line:
                        istep += 1
                    if "Total Energy =" in line:
                        word = line.split()
                        toten[iz, scan_path[istep][0], scan_path[istep][1]] = float(word[3]) * rydberg

    # write toten.dat
    with open("toten.dat", "w") as f2:
        f2.write("#   ix    iy    iz      toten(eV)\n")
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    f2.write(f"{ix:6d}{iy:6d}{iz:6d}{toten[iz][iy][ix]:24.12f}\n")

# ==================== caclulate forces for tilt correction ====================
if use_tilt:
    # interpolation from a 2d grid
    toten_1d = []
    for iz in range(nz):
        toten_1d.append(scipy.interpolate.UnivariateSpline(y_grid, toten[iz, :, 0], s=0))  # cubic spline

    # We use an iterative method to solve D=F(D)/k. For conventional D=F/k, set Niter=1 and alpha=1
    y_init = y_grid[np.newaxis, :]
    y_new = np.tile(y_init, (nz, 1))
    for iz in range(nz):
        for iiter in range(niter):
            fy = (y_grid - y_new[iz, :]) * k_spring
            fy += (toten_1d[iz](y_new[iz, :] - h*0.5) - toten_1d[iz](y_new[iz, :] + h*0.5)) / h
            y_incr = fy / k_spring
            y_new[iz, :] += y_incr * alpha
            if np.max(np.abs(y_incr)) < 1.e-4:
                if verbose:
                    print(f" iz: {iz:0d}, tilt correction converge at iiter = {iiter:0d}")
                break
            if iiter == niter-1:
                print(f"Warning: iz: {iz:0d}, tilt correction does not converge in niter = {niter:0d}")

# ==================== calculate kts ====================
if plot_force:
    # in the plot of force, the variable kts is the force.
    if use_tilt:
        kts = np.zeros((nz-2, ny))
        for iz in range(1, nz-1):
            kts[iz-1, :] = (toten_1d[iz-1](y_new[iz, :]) - toten_1d[iz+1](y_new[iz, :])) / (2*z_spacing)
    else:
        kts = (toten[0:nz-2, :, 0] - toten[2:nz, :, 0]) / (2*z_spacing)

else:
    if use_tilt:
        kts = np.zeros((nz-2, ny))
        for iz in range(1, nz-1):
            kts[iz-1, :] = (toten_1d[iz-1](y_new[iz, :]) - 2*toten_1d[iz](y_new[iz, :])
                            + toten_1d[iz+1](y_new[iz, :])) / z_spacing**2
    else:
        kts = (toten[0:nz-2, :, 0] - 2*toten[1:nz-1, :, 0] + toten[2:nz, :, 0]) / z_spacing**2

if use_bohr:
    if plot_force:
        # convert k_ts from eV/A to Ha/a0
        kts = kts * bohr/Ha
    else:
        # convert k_ts from eV/A^2 to Ha/a0^2
        kts = kts * bohr**2/Ha

# ==================== calculate maximums ====================
if show_max:
    interp_func = scipy.interpolate.RectBivariateSpline(z_grid[1:nz-1], y_grid, kts, kx=3, ky=3)
    z_max_list = np.linspace(z_grid[1], z_grid[nz-2], nz-2)
    y_max_list = []
    for z in z_max_list:
        def y_grid_1d(y): return -interp_func(z, y)[0]  # Negative for maximizing
        z_max_result = scipy.optimize.minimize_scalar(y_grid_1d, bounds=(y_range[0], y_range[1]), method='bounded')
        y_max_list.append(z_max_result.x[0])
    print(" iz    z(A)    y(A)")
    for iz in range(nz-2):
        print(f"{iz+1:3d}{z_max_list[iz]:8.4f}{y_max_list[iz]:8.4f}")

# ==================== construct the atomic structure ====================
palette = load_palette()
if use_bohr:
    funit = bohr
else:
    funit = 1.0
if show_atom:
    color_dict = load_atom_color()
    poscar1 = Poscar()
    poscar1.read_parsec(filename="sample.parsec_st.dat")
    # poscar1 is modified here. If a feature needs to keep poscar1, use copy.
    if poscar1.Ndim == 2:
        poscar1.supercell([2, 2, 1])
        poscar1.move([-0.5, -0.5, 0], lbox=False)
    apc = poscar1.cartesian()
    atom = poscar1.atom_list()

    zmax = np.max(apc[:, 2])

    atom_y = []
    atom_color = []
    atom_name = []
    for ia in range(len(apc)):
        if apc[ia, 2] > zmax - 1.0 and abs(apc[ia, 0]) < 1.e-4:
            atom_y.append(apc[ia, 1]/funit)
            atom_color.append(palette[color_dict[atom[ia]]])
            atom_name.append(atom[ia])

# ==================== creating the plot ====================
mpl.rcParams["font.sans-serif"].insert(0, "Noto Sans")
mpl.rcParams.update({'font.size': 14})
mpl.rcParams.update({'mathtext.default': 'regular'})

fig0 = plt.figure(figsize=(5, 3.75))
gs0 = fig0.add_gridspec(1, 2, wspace=0.02, hspace=0.00, left=0.14, right=0.80,
                        top=0.95, bottom=0.15, width_ratios=[0.6, 0.04])
[ax0, ax1] = gs0.subplots()

im_extent = [y_range[0]-y_spacing*0.5, y_range[1]+y_spacing*0.5, z_grid[1]-z_spacing*0.5, z_grid[nz-2]+z_spacing*0.5]
for ic in range(len(im_extent)):
    im_extent[ic] = im_extent[ic]/funit
if plot_force:
    norm = mpl.colors.TwoSlopeNorm(vmin=-0.21, vcenter=0, vmax=0.42)  # may depend on the system: use kts.max()
    cmap = "seismic"
    im = ax0.imshow(kts, interpolation='spline36', cmap=cmap, origin="lower", norm=norm, extent=im_extent, zorder=1)
else:
    im = ax0.imshow(kts, interpolation='spline36', cmap="rainbow", origin="lower", extent=im_extent, zorder=1)

ax0.set_aspect((y_range[1]-y_range[0])/(z_range[1]-z_range[0]-z_spacing*2))

if show_atom:
    for ia in range(len(atom_y)):
        ax0.axvline(x=atom_y[ia], color=atom_color[ia], linestyle="dashdot", linewidth=1)
if show_max:
    ax0.plot(y_max_list, z_max_list, color=palette["red"], linestyle="dotted", linewidth=1)

ax0.set_xlim([y_range[0]/funit, y_range[1]/funit])
ax0.set_ylim([z_grid[1]/funit, z_grid[nz-2]/funit])
if use_bohr:
    ax0.set_xlabel(r"$\mathit{y}\ (Bohr)$", color=palette["black"])
    ax0.set_ylabel(r"$\mathit{z}\ (Bohr)$", color=palette["black"])
else:
    ax0.set_xlabel(r"$\mathit{y}\ (Å)$", color=palette["black"])
    ax0.set_ylabel(r"$\mathit{z}\ (Å)$", color=palette["black"])

if plot_force:
    cb = fig0.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax1, orientation='vertical', extend='both')
else:
    cb = fig0.colorbar(im, cax=ax1, orientation='vertical')

cb.outline.set_linewidth(1)
cb.outline.set_color(palette["black"])
if plot_force:
    if use_bohr:
        ax1.set_ylabel(r"$\mathit{F}\ (a.u.)$", color=palette["black"])
    else:
        ax1.set_ylabel(r"$\mathit{F}\ (eV/Å)$", color=palette["black"])
else:
    if use_bohr:
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
if plot_force:
    filename += "_force"
if use_tilt:
    filename += "_tilt"
if show_atom:
    filename += "_atom"
if show_max:
    filename += "_max"
if use_bohr:
    filename += "_bohr"
filename += ".png"
fig0.savefig(filename, dpi=1200)

# ==================== vector map for tilt corrections ====================
if use_tilt:
    fig1 = plt.figure(figsize=(5, 3.75))
    gs1 = fig1.add_gridspec(1, 1, left=0.14, right=0.74, top=0.95, bottom=0.15)
    ax2 = gs1.subplots()

    # q = ax2.quiver(y_grid/funit, z_grid/funit, (y_new-y_init)/funit, np.zeros((nz,ny)), angles='xy',
    #               scale_units='xy', scale=1, color=palette["darkblue"], linewidth=1, zorder=3)
    for iy in range(ny):
        ax2.plot(y_new[:, iy], z_grid, color=palette["darkblue"], linewidth=1, zorder=3)

    if show_atom:
        for ia in range(len(atom_y)):
            ax2.axvline(x=atom_y[ia], color=atom_color[ia], linestyle="dashdot", linewidth=1)
    ax2.set_xlim([y_range[0]/funit, y_range[1]/funit])
    ax2.set_ylim([z_range[0]/funit, z_range[1]/funit])
    if use_bohr:
        ax2.set_xlabel(r"$\mathit{y}\ (Bohr)$", color=palette["black"])
        ax2.set_ylabel(r"$\mathit{z}\ (Bohr)$", color=palette["black"])
    else:
        ax2.set_xlabel(r"$\mathit{y}\ (Å)$", color=palette["black"])
        ax2.set_ylabel(r"$\mathit{z}\ (Å)$", color=palette["black"])

    ax2.tick_params(axis="x", bottom=True, right=False, direction="in",
                    color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
    ax2.tick_params(axis="y", left=True, right=False, direction="in",
                    color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0)
    for edge in ["bottom", "top", "left", "right"]:
        ax2.spines[edge].set_color(palette["black"])
        ax2.spines[edge].set_linewidth(1)
        ax2.spines[edge].set_zorder(4)

    filename = "tilt"
    if show_atom:
        filename += "_atom"
    if use_bohr:
        filename += "_bohr"
    filename += ".png"
    fig1.savefig(filename, dpi=1200)
