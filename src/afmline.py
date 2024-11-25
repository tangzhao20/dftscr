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

f1 = open("afm.in", "r")
line = f1.readlines()
f1.close()
for l in line:
    word = l.split()
    if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
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
nx = 0
xlist = []
while (x < x_range[1]+1e-6):
    xlist.append(x)
    x += x_spacing
    nx += 1
y = y_range[0]
ny = 0
ylist = []
while (y < y_range[1]+1e-6):
    ylist.append(y)
    y += y_spacing
    ny += 1
z = z_range[0]
nz = 0
zlist = []
while (z < z_range[1]+1e-6):
    zlist.append(z)
    z += z_spacing
    nz += 1

f3 = open("steps.dat", "r")
line = f3.readlines()
f3.close()
steplist = []
for l in line:
    word = l.split()
    if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
        continue
    steplist.append(int(word[0]))

lxincrease = True
k = 0
ip = 0
# movelist here is the index (ix,iy) of each point
movelist = []
for iy in range(ny):
    if lxincrease:
        for ix in range(nx):
            if k == 0:
                movelist.append([])
            movelist[-1].append([ix, iy])
            k += 1
            if k == steplist[ip]:
                ip += 1
                k = 0
    else:
        for ix in range(nx-1, -1, -1):
            if k == 0:
                movelist.append([])
            movelist[-1].append([ix, iy])
            k += 1
            if k == steplist[ip]:
                ip += 1
                k = 0
    lxincrease = not lxincrease

# ==================== calculate or read toten ====================

files = os.listdir()
# if the toten.dat exist, read it, if not, read from calculations outputs
if "toten.dat" in files:
    f5 = open("toten.dat", "r")
    line = f5.readlines()
    f5.close()
    il = 1
    toten = []
    for iz in range(nz):
        toten1 = []
        for iy in range(ny):
            toten0 = []
            for ix in range(nx):
                toten0.append(float(line[il].split()[3]))
                il += 1
            toten1.append(toten0)
        toten.append(toten1)

else:
    # initialize toten matrix
    # toten[nz][ny][nx] in Ry
    toten = []
    for iz in range(nz):
        toten0 = []
        for iy in range(ny):
            toten0.append([0.0]*nx)
        toten.append(toten0)

    for iz in range(nz):
        for ip in range(parallel):
            istep = -1
            filename = "seq_"+str(iz+1)+"_"+str(ip+1)+"/parsec.out"
            f1 = open(filename, "r")
            line = f1.readlines()
            f1.close()

            for l in line:
                word = l.split()
                if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                    continue
                if len(word) >= 2 and word[0] == "Starting" and word[1] == "SCF...":
                    istep += 1
                if len(word) >= 5 and word[0] == "Total" and word[1] == "Energy" and word[2] == "=":
                    toten[iz][movelist[ip][istep][1]][movelist[ip][istep]
                                                      [0]] = float(word[3])*rydberg  # convert Ry to eV

    # write the toten file
    f5 = open("toten.dat", "w")
    f5.write("#ix iy iz toten(eV)\n")
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                f5.write(str(ix)+" "+str(iy)+" "+str(iz)+" "+str(toten[iz][iy][ix])+"\n")
    f5.close()

# ==================== caclulate kts ====================

kts = []
ix = 0
for iz in range(1, nz - 1):
    kts0 = []
    for iy in range(ny):
        kts1 = (toten[iz-1][iy][0] - 2*toten[iz][iy][0] + toten[iz+1][iy][ix]) / z_spacing**2
        if lbohr:
            # convert k_ts from eV/A^2 to Ha/a0^2
            kts1 = kts1 * bohr**2 / Ha
        kts0.append(kts1)
    kts.append(kts0)

# ==================== caclulate maximums ====================

if lmax:
    interp_func = scipy.interpolate.RectBivariateSpline(zlist[1:nz-1], ylist, kts, kx=3, ky=3)
    z_max_list = np.linspace(z_range[0], z_range[1], 500)
    y_max_list = []
    for z in z_max_list:
        def ylist_1d(y): return -interp_func(z, y)[0]  # Negative for maximizing
        z_max_result = scipy.optimize.minimize_scalar(ylist_1d, bounds=(y_range[0], y_range[1]), method='bounded')
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

fig = plt.figure(figsize=(5, 3.75))
gs0 = fig.add_gridspec(1, 2, wspace=0.02, hspace=0.00, left=0.14, right=0.80,
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

cb = fig.colorbar(im, cax=ax1, orientation='vertical')
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
fig.savefig(filename, dpi=1200)
