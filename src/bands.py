#!/usr/bin/env python3

# bands.py package (E1) (E2)
# Making the band structure plot

# if E2 exists, the energy range is (E1,E2)
# if only E1 exists, the energy range is (-E1,E1)
# if neither exists, the energy range is (-5 eV, 5 eV)

import os
import sys
from classes import *
from load_data import load_package_name

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
        eigenval1.eig -= eigenval1.vbm
        eigenval1.writegap(kpoints1)
    else:
        doscar1 = Doscar()
        doscar1.read_vasp()
        eigenval1.eig -= doscar1.ef

elif package in package_name["qe"]+package_name["qeproj"]:
    # Input: *.xml, kpath.in
    # No need to run bands.x

    poscar1.read_qe_xml()
    rlc = poscar1.rlc()

    eigenval1.read_qe()

    kpoints1.read_kpathin()

    if eigenval1.is_semic == True:
        eigenval1.eig -= eigenval1.vbm
        eigenval1.writegap(kpoints1)
    else:
        doscar1 = Doscar()
        doscar1.read_qe_xml()
        eigenval1.eig -= doscar1.ef

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

    eigenval1.read_wannier90(Nb_pad=pad)

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
        eigenval2.read_qe("../bands/" + filename)

        eigenval2.eig -= eigenval2.vbm
        if eigenval2.is_semic == True:
            eigenval1.is_semic = True

    eigenval1.occ = np.zeros((eigenval1.Ns, eigenval1.Nk, eigenval1.Nb))
    eigenval1.occ[:, :, :eigenval2.Nvb[0] - pad] = 1.0

    eigenval1.calculate_gap()
    if eigenval1.is_semic == True:
        eigenval1.eig -= eigenval1.vbm
        eigenval1.writegap(kpoints1)
    else:
        print("Metal band structure are not shifted")

elif package in package_name["parsec"]:
    # Input: bands.dat, parsec.in, kpath.in

    poscar1.read_parsec()
    rlc = poscar1.rlc()

    eigenval1.read_parsec(lc=poscar1.lc)

    kpoints1.read_kpathin()

    if eigenval1.is_semic == True:
        eigenval1.eig -= eigenval1.vbm
        eigenval1.writegap(kpoints1)
    else:
        print("Metal band structure are not shifted")


else:
    print("Package \""+package+"\" is not supported yet.")
    print("python bands.py package (Emin) (Emax)")
    sys.exit()

is_proj = False
if package in package_name["vaspproj"]+package_name["qeproj"]:
    is_proj = True

    if package in package_name["vaspproj"]:
        # Input: PROCAR
        eigenval1.read_vasp_procar()
    elif package in package_name["qeproj"]:
        # Input: projwfc.out
        eigenval1.read_qe_projwfc()

    atom_list = sys.argv[2]
    del sys.argv[2]
    atom_flag = poscar1.read_atom_list(atom_list)
    orb_list = sys.argv[2]
    del sys.argv[2]
    orb_flag = eigenval1.read_orb_list(orb_list)

# The following plotting section should be generalized, not code-specific

x_ticks = kpoints1.x_ticks_out(rlc)
x_labels = kpoints1.x_labels_out()

y_range = [-5.0, 5.0]
if len(sys.argv) >= 4:
    y_range[0] = min(float(sys.argv[2]), float(sys.argv[3]))
    y_range[1] = max(float(sys.argv[2]), float(sys.argv[3]))
elif len(sys.argv) == 3:
    y_range[0] = -float(sys.argv[2])
    y_range[1] = float(sys.argv[2])

bands = BandsPlot(x_ticks, x_labels, y_range)
bands.plot_bands(eigenval1, kpoints1, rlc)

output_name = "bs.png"

# projection plot
if is_proj:
    if eigenval1.Ns == 1:
        proj_color = ["orange"]
    else:
        proj_color = ["blue", "beige"]
    dot_size = 50.0
    proj_plot_size = eigenval1.plot_proj(atom_flag, orb_flag) * dot_size
    for ispin in range(eigenval1.Ns):
        bands.add_scatter(bands.x, bands.energy[ispin], proj_plot_size[ispin, :, :],
                          color=bands.palette[proj_color[ispin]], zorder=2.5-ispin)

    output_name = "proj_"+atom_list+"_"+orb_list+".png"

    f2 = open("proj_"+atom_list+"_"+orb_list+".dat", "w")
    f2.write("#ispin x_kpath energy(eV) size\n")
    for ispin in range(eigenval1.Ns):
        for ib in range(eigenval1.Nb):
            ik0 = 0
            for ip in range(bands.Np):
                for ik in range(len(bands.x[ip])):
                    f2.write(f"{ispin}  {bands.x[ip][ik]}  {bands.energy[ispin][ib][ik0]}  " +
                             f"{proj_plot_size[ispin, ib, ik0]}\n")
                    ik0 += 1
            f2.write("\n")
    f2.close()


# Second band structure plot (for wannier)
if fsecond:
    x2 = eigenval2.eig_x(kp=kpoints1, rlc=rlc)
    energy2 = eigenval2.eig.swapaxes(1, 2)  # (Ns, Nk, Nb) -> (Ns, Nb, Nk)
    bands.add_plot(x2, energy2[0], color=bands.palette["orange"], label="Wannier", zorder=2)
    bands.write_bands(x=x2, energy=energy2, filename="eigenval_wan.dat")

bands.fig.savefig(output_name, dpi=1200)
