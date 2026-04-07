#!/usr/bin/env python3

"""
This script works with afm.sh together to prepare the inputs for AFM simulation

inputs: afm.in, tip.xyz, sample.parsec_st.dat, *_POTRE.DAT

example:
    python afm.py --vasp
"""

import argparse
import sys

import numpy as np

from classes import Poscar
from load_data import load_constant, load_atom_index


bohr = load_constant("bohr")

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument("--vasp", "-v", action="store_true", help="export the structure to afm.vasp for visualization")
parser.add_argument("--verbose", action="store_true", help="print verbose output")

args = parser.parse_args()

export_vasp = args.vasp
verbose = args.verbose

# ================================================================
# read the input file
z_range = [5.7, 6.3]
x_spacing = 0.6
y_spacing = 0.6
z_spacing = 0.3
boundary = -1e0
vacuum = 7.5
parallel = 1
use_fdet = False
use_spin = False

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

            case ["boundary", val, *rest]:
                if rest and rest[0].lower().startswith("a"):
                    boundary = float(val) / bohr
                else:
                    boundary = float(val)
            case ["vacuum", val]: vacuum = float(val)

            case ["parallel", val]: parallel = int(val)

            case ["fdet", *rest]:
                use_fdet = not (rest and rest[0].lower() in ["false", ".false."])
            case ["spin", *rest]:
                use_spin = not (rest and rest[0].lower() in ["false", ".false.", "1", "0"])

            case [("contrast_range" | "k_spring" | "niter" | "alpha" | "h"), *_]:
                pass

            case _:
                print(f"Warning: check your input line: {' '.join(word)}")

if verbose:
    print(f"x_range: {x_range[0]:f} ~ {x_range[1]:f} Bohr")
    print(f"y_range: {y_range[0]:f} ~ {y_range[1]:f} Bohr")
    print(f"z_range: {z_range[0]:f} ~ {z_range[1]:f} Bohr")
    print(f"x_spacing: {x_spacing:f} Bohr")
    print(f"y_spacing: {y_spacing:f} Bohr")
    print(f"z_spacing: {z_spacing:f} Bohr")
    if use_fdet:
        print("FDET is enabled")
    else:
        print("FDET is not enabled")
    if use_spin:
        print("Nspin is set to 2")
    else:
        print("Nspin is set to 1")

x_grid = np.arange(x_range[0], x_range[1]+1e-6, x_spacing)
nx = len(x_grid)
x_range[1] = x_grid[-1]

y_grid = np.arange(y_range[0], y_range[1]+1e-6, y_spacing)
ny = len(y_grid)
y_range[1] = y_grid[-1]

z_grid = np.arange(z_range[0], z_range[1]+1e-6, z_spacing)
nz = len(z_grid)
z_range[1] = z_grid[-1]

# ================================================================
# Read the structure from tip.xyz and sample.parsec_st.dat
poscar1 = Poscar()
poscar1.read_xyz("tip.xyz")
poscar2 = Poscar()
poscar2.read_parsec("sample.parsec_st.dat")

if poscar2.Ndim != 0 and poscar2.Ndim != 2:
    print("Error: Only Ndim = 0 or 2 are supported.")
    sys.exit()

# ================================================================
# Read the core charge from pp files *_POTRE.DAT
zion1 = []
for i in range(poscar1.Ntype):
    pp_name = f"{poscar1.atomtype[i]}_POTRE.DAT"
    with open(pp_name, "r") as f5:
        lines = f5.readlines()
        zion1.append(int(float(lines[3].split()[5]) + 0.5))
zion2 = []
for i in range(poscar2.Ntype):
    pp_name = f"{poscar2.atomtype[i]}_POTRE.DAT"
    with open(pp_name, "r") as f5:
        lines = f5.readlines()
        zion2.append(int(float(lines[3].split()[5]) + 0.5))

# Calculate the maximum of Nb for in advance.
atomdir = load_atom_index()
z_full_list = [2, 10, 18, 36, 54, 86, 118]
Nb1_max = 0
for i in range(poscar1.Ntype):
    # Nb_atom=(z_full-z_atom+z_pp)/2
    z_atom = atomdir[poscar1.atomtype[i]]
    for j in range(7):
        if z_full_list[j] >= z_atom:
            z_full = z_full_list[j]
            break
    z_pp = zion1[i]
    Nb1_max += (z_full-z_atom+z_pp)*poscar1.Naint[i]
Nb1_max = Nb1_max//2
Nb2_max = 0
for i in range(poscar2.Ntype):
    # Nb_atom=(z_full-z_atom+z_pp)/2
    z_atom = atomdir[poscar2.atomtype[i]]
    for j in range(7):
        if z_full_list[j] >= z_atom:
            z_full = z_full_list[j]
            break
    z_pp = zion2[i]
    Nb2_max += (z_full-z_atom+z_pp)*poscar2.Naint[i]
Nb2_max = Nb2_max//2
# Calculate the prefered Nb of sample.
Nb2 = 0
for i in range(poscar2.Ntype):
    Nb2 += zion2[i]*poscar2.Naint[i]
# Modify the following line if needed.
Nb2 = int(Nb2/2*1.1+10)
Nb2 = min(Nb2, Nb2_max)

# calculate the apc, this should be done before the loop
apc1 = poscar1.cartesian(factor=1.0/bohr)
apc2 = poscar2.cartesian(factor=1.0/bohr)

# move the tip atom to (0,0,0)
itip = np.argmin(apc1[:, 2])
apc1 = apc1 - apc1[itip, :]

# move the sample top to z=0
zmax = np.max(apc2[:, 2])

# instead of the very top atom, we select the medium of atoms with in 1 A under the top
z_toplist = np.sort(apc2[apc2[:, 2] > (zmax - 1.0 / bohr), 2])
imid = len(z_toplist) // 2
zmid = (z_toplist[imid] + z_toplist[~imid]) / 2
apc2[:, 2] -= zmid

# ================================================================
# calculate the step numbers. Need to minus 1 in the parsec.in file
steps_per_task = [(nx * ny) // parallel] * parallel
for ip in range((nx * ny) % parallel):
    steps_per_task[ip] += 1
with open("steps.dat", "w") as f3:
    for steps in steps_per_task:
        f3.write(f"{steps}\n")
task_paths = []

is_x_forward = True
k = 0
ip = 0
for iy in range(ny):
    if is_x_forward:
        for ix in range(nx):
            if k == 0:
                task_paths.append([])
            task_paths[-1].append([x_grid[ix], y_grid[iy]])
            k += 1
            if k == steps_per_task[ip]:
                ip += 1
                k = 0
    else:
        for ix in range(nx-1, -1, -1):
            if k == 0:
                task_paths.append([])
            task_paths[-1].append([x_grid[ix], y_grid[iy]])
            k += 1
            if k == steps_per_task[ip]:
                ip += 1
                k = 0
    is_x_forward = not is_x_forward

# ================================================================
# If the boundary is not set, calculate it
if boundary < 0:
    if poscar2.Ndim == 0:
        for ix0 in range(2):
            for ix1 in range(2):
                offset = np.array([x_range[ix0], y_range[ix1], z_range[1]])
                boundary = max(boundary, np.linalg.norm(apc1+offset, axis=1).max())
        boundary = max(boundary, np.linalg.norm(apc2, axis=1).max())
        boundary += vacuum

    elif poscar2.Ndim == 2:
        boundary_max = apc1[:, 2].max() + z_range[1]
        boundary_min = apc2[:, 2].min()
        # Move the atoms to centralize the sample-tip structure
        z_move = 0.5 * (boundary_max + boundary_min)
        apc1[:, 2] -= z_move
        apc2[:, 2] -= z_move
        boundary = 0.5 * (boundary_max - boundary_min) + vacuum

else:  # Calculate the z_move, make the bottom at 10 bohr from boundary
    if poscar2.Ndim == 2:
        boundary_min = apc2[:, 2].min()
        z_move = boundary_min - (-boundary + vacuum)
        apc1[:, 2] -= z_move
        apc2[:, 2] -= z_move


# if FDET, we do a calculation for a sample potential file
if use_fdet:
    with open("parsec_st_spot.dat", "w") as f2:
        f2.write("#---------output from afm.py----------\n")
        if use_spin:
            f2.write("spin_polarization .true.\n")
        f2.write(f"states_num {Nb2}\n\n")
        if poscar2.Ndim == 0:
            pass
        if poscar2.Ndim == 2:
            f2.write("boundary_conditions slab\n")
            f2.write("begin cell_shape\n")
            for ix0 in range(2):
                for ix1 in range(3):
                    f2.write(f"{poscar2.lc[ix0, ix1]/bohr:18.12f}")
                f2.write("\n")
            f2.write("end cell_shape\n\n")

            f2.write("kpoint_method mp\n\n")
            f2.write("begin monkhorst_pack_grid\n")
            k_grid = poscar2.find_k_grid()
            for ix in range(2):
                f2.write(f"  {k_grid[ix]:d}")
            f2.write("\nend monkhorst_pack_grid\n\n")
            k_grid_shift = poscar2.find_k_grid_shift()
            f2.write("begin monkhorst_pack_shift\n")
            for ix in range(2):
                f2.write(f"  {k_grid_shift[ix]:f}")
            f2.write("\nend monkhorst_pack_shift\n\n")

        f2.write(f"boundary_sphere_radius {boundary:.12g}\n\n")
        f2.write(f"atom_types_num {poscar2.Ntype}\n")
        f2.write("coordinate_unit cartesian_bohr\n\n")

        f2.write("#------------ begin sample -----------\n")
        k = 0
        for i in range(poscar2.Ntype):
            f2.write(f"atom_type {poscar2.atomtype[i]}\n")
            f2.write("local_component s\n")
            f2.write("begin atom_coord\n")
            for ia in range(poscar2.Naint[i]):
                f2.write(f"{apc2[k, 0]:18.12f}{apc2[k, 1]:18.12f}{apc2[k, 2]:18.12f}\n")
                k = k+1
            f2.write("end atom_coord\n\n")
        f2.write("#------------- end sample ------------\n\n")

# main loop
for iz in range(nz):
    for ip in range(parallel):
        # ================================================================
        # Write the structure to parsec_st_iz_ip.dat file
        with open(f"parsec_st_{iz+1}_{ip+1}.dat", "w") as f2:
            f2.write("#---------output from afm.py----------\n\n")
            f2.write("minimization manual\n\n")
            if use_spin and not use_fdet:
                f2.write("spin_polarization .true.\n")
            if use_fdet:
                f2.write(f"states_num {Nb1_max}\n\n")
            else:
                f2.write(f"states_num {Nb1_max+Nb2}\n\n")
            if poscar2.Ndim == 0:
                pass
            if poscar2.Ndim == 2:
                f2.write("boundary_conditions slab\n")
                f2.write("begin cell_shape\n")
                for ix0 in range(2):
                    for ix1 in range(3):
                        f2.write(f"{poscar2.lc[ix0, ix1]/bohr:18.12f}")
                    f2.write("\n")
                f2.write("end cell_shape\n\n")

                f2.write("kpoint_method mp\n\n")
                f2.write("begin monkhorst_pack_grid\n")
                k_grid = poscar2.find_k_grid()
                for ix in range(2):
                    f2.write(f"  {k_grid[ix]:d}")
                f2.write("\nend monkhorst_pack_grid\n\n")
                k_grid_shift = poscar2.find_k_grid_shift()
                f2.write("begin monkhorst_pack_shift\n")
                for ix in range(2):
                    f2.write(f"  {k_grid_shift[ix]:f}")
                f2.write("\nend monkhorst_pack_shift\n\n")

            f2.write(f"boundary_sphere_radius {boundary:.12g}\n\n")
            if use_fdet:
                f2.write(f"atom_types_num {poscar1.Ntype}\n")
            else:
                f2.write(f"atom_types_num {poscar1.Ntype+poscar2.Ntype}\n")
            f2.write("coordinate_unit cartesian_bohr\n\n")

            f2.write("#------------- begin tip -------------\n")
            k = 0
            for i in range(poscar1.Ntype):
                f2.write(f"atom_type {poscar1.atomtype[i]}\n")
                # Assume the local_component is s here. A dictionary of {element: local_component} should be added
                f2.write("local_component: s\n")
                f2.write("begin atom_coord\n")
                for ia in range(poscar1.Naint[i]):
                    f2.write(f"{apc1[k, 0]+task_paths[ip][0][0]:18.12f}{apc1[k, 1]+task_paths[ip][0][1]:18.12f}" +
                             f"{apc1[k, 2]+z_grid[iz]:18.12f}\n")
                    k = k+1
                f2.write("end atom_coord\n\n")
            f2.write("#-------------- end tip --------------\n\n")

            f2.write("#------------ begin sample -----------\n")
            if use_fdet:
                f2.write("add_point_charges .TRUE.\n")
                f2.write(f"point_typ_num {poscar2.Ntype}\n\n")

            k = 0
            for i in range(poscar2.Ntype):
                if use_fdet:
                    f2.write(f"pt_chg: {zion2[i]}\n")
                    f2.write("begin point_coord\n")
                else:
                    f2.write(f"atom_type: {poscar2.atomtype[i]}\n")
                    f2.write("local_component: s\n")
                    f2.write("begin atom_coord\n")

                for ia in range(poscar2.Naint[i]):
                    f2.write(f"{apc2[k, 0]:18.12f}{apc2[k, 1]:18.12f}{apc2[k, 2]:18.12f}\n")
                    k = k+1

                if use_fdet:
                    f2.write("end point_coord\n\n")
                else:
                    f2.write("end atom_coord\n\n")
            f2.write("#------------- end sample ------------\n\n")

        # write the manual_*.dat file
        with open(f"manual_{iz+1}_{ip+1}.dat", "w") as f4:
            for istep in range(1, steps_per_task[ip]):
                dx = task_paths[ip][istep][0]
                dy = task_paths[ip][istep][1]
                dz = z_grid[iz]
                for ia in range(poscar1.Natom):
                    f4.write(f"{apc1[ia, 0] + dx:18.12f}{apc1[ia, 1] + dy:18.12f}{apc1[ia, 2] + dz:18.12f}\n")
                if not use_fdet:
                    for ia in range(poscar2.Natom):
                        f4.write(f"{apc2[ia, 0]:18.12f}{apc2[ia, 1]:18.12f}{apc2[ia, 2]:18.12f}\n")
                f4.write("\n")

# convert the structure to vasp format
if export_vasp:
    poscar3 = Poscar()
    if use_fdet:
        poscar3.read_parsec("parsec_st_spot.dat")
    else:
        poscar3.read_parsec("parsec_st_1_1.dat")
    pi = load_constant("pi")
    rlc3 = poscar3.rlc()
    if poscar2.Ndim == 0:
        shift = [0.5, 0.5, 0.5]
    elif poscar2.Ndim == 2:
        shift = [0.0, 0.0, 0.5]
    shift = np.array(shift)
    offset = np.array([sum(x_range)/2, sum(y_range)/2, sum(z_range)/2])
    ap1 = (apc1 + offset) * bohr @ rlc3 / (2*pi) + shift
    if use_fdet:
        ia = 0
        for itype in range(poscar1.Ntype):
            for iaint in range(poscar1.Naint[itype]):
                poscar3.add_atom(itype=itype, ap=ap1[ia, :], new_type=poscar1.atomtype[itype])
                ia += 1
    else:
        poscar3.ap[:poscar1.Natom, :] = ap1
    poscar3.write_vasp("afm.vasp")
