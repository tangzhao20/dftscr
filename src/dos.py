#!/usr/bin/env python3

# plot DOS

# python dos.py (v) package (E1) (E2)

# VASP input: DOSCAR
# QE input: *.dos *.xml
# QE projection input: *.dos *.xml *.pdos_atm#*(*)_wfc#*(*)

import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from classes import Doscar, Poscar, Procar
from load_data import load_package_name, load_palette

fvertical = False
filename_out = "dos.png"

for iw in range(1, len(sys.argv)):
    if sys.argv[iw] in ["v", "vertical"]:
        fvertical = True
        filename_out = "dos_v.png"
        del sys.argv[iw]
        break
if len(sys.argv) <= 1:
    print("python dos.py package (v) (Emin) (Emax)")
    sys.exit()
package = sys.argv[1]

package_name = load_package_name()

doscar0 = Doscar()
if package in package_name["vasp"]+package_name["vaspproj"]:
    doscar0.read_vasp()
elif package in package_name["qe"]+package_name["qeproj"]:
    # find a *.dos file
    files = os.listdir(".")
    fxml = False
    for f in files:
        if f.endswith('.dos'):
            filename = f
        if f.endswith('.xml'):  # sometimes EFermi in .dos file could be wrong??
            fxml = True          # read from .xml file instead.
            xmlname = f
    doscar0.read_qe(filename)
    if fxml:
        print("EFermi read from .xml file")
        doscar0.read_xml(xmlname)
else:
    print("Package \""+package+"\" is not supported yet.")
    print("python dos.py (v) package (Emin) (Emax)")
    sys.exit()

is_proj = False
input_shift = 0
if package in package_name["vaspproj"]+package_name["qeproj"]:
    is_proj = True
    input_shift = 2

if len(sys.argv) >= 4 + input_shift:
    xmax = float(sys.argv[3 + input_shift])
    xmin = float(sys.argv[2 + input_shift])
elif len(sys.argv) == 3 + input_shift:
    xmax = float(sys.argv[2 + input_shift])
    xmin = -float(sys.argv[2 + input_shift])
else:
    xmax = 5.0
    xmin = -5.0

if is_proj:
    poscar0 = Poscar()
    procar0 = Procar()
    if package in package_name["vaspproj"]:
        sigma = 0.05
        poscar0.read_vasp()
        procar0.read_vasp()
        doscar0.Nepdos = 1000
        doscar0.energy_pdos = np.linspace(xmin+doscar0.ef, xmax+doscar0.ef, doscar0.Nepdos)
        doscar0.pdos = procar0.calculate_pdos(doscar0.energy_pdos, sigma)
        doscar0.has_pdos = True

    if package in package_name["qeproj"]:
        # qe support is broken in this commit
        # procar and poscar should be initialized here
        doscar0.readpdos_qe()

    atom_list = sys.argv[2]
    del sys.argv[2]
    atom_flag = poscar0.read_atom_list(atom_list)

    orb_list = sys.argv[2]
    del sys.argv[2]
    orb_flag = procar0.read_orb_list(orb_list)

    pdos = doscar0.plot(atom_flag, orb_flag)

# The following plotting section should be generalized, not code-specific

doscar0.energyshift(doscar0.ef)

dosmax = 0.0
for ie in range(doscar0.Nedos):
    if doscar0.energy[ie] < xmin:
        continue
    if doscar0.energy[ie] > xmax:
        break
    if doscar0.Ns == 1:
        dosmax = max(dosmax, doscar0.dos[0][ie])
    else:
        dosmax = max(dosmax, doscar0.dos[0][ie], doscar0.dos[1][ie])

palette = load_palette()
mpl.rcParams["font.sans-serif"].insert(0, "Noto Sans")
mpl.rcParams.update({'font.size': 14})

if doscar0.Ns == 1:
    color_dos = ["darkblue"]
    color_pdos = ["orange"]
else:
    color_dos = ["darkblue", "orange"]
    color_pdos = ["blue", "beige"]

if fvertical == False:
    fig = plt.figure(figsize=(5, 3.75))
    gs0 = fig.add_gridspec(doscar0.Ns, 1, wspace=0.0, hspace=0.00, left=0.14, right=0.96, top=0.97, bottom=0.12)
    ax = []
    for ispin in range(doscar0.Ns):
        ax.append(fig.add_subplot(gs0[ispin]))

    for ispin in range(doscar0.Ns):
        ax[ispin].axvline(linewidth=1, color=palette["gray"], zorder=0)
        ax[ispin].plot(doscar0.energy, doscar0.dos[ispin], color=palette[color_dos[ispin]], linewidth=1, zorder=3)

    if is_proj:
        for ispin in range(doscar0.Ns):
            ax[ispin].plot(doscar0.energy_pdos, pdos[ispin, :], color=palette[color_pdos[ispin]],
                           linewidth=1, zorder=3.5)
        filename_out = "dos_"+atom_list+"_"+orb_list+".png"

    ax[0].set_ylim([0, dosmax*1.1])
    if doscar0.Ns == 2:
        ax[1].set_ylim([dosmax*1.1, 0])
    for ispin in range(doscar0.Ns):
        ax[ispin].set_xlim([xmin, xmax])
        ax[ispin].tick_params(axis="x", bottom=True, top=True, direction="in",
                              color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0, pad=4)
        ax[ispin].tick_params(axis="y", left=True, right=True, direction="in",
                              color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0, pad=4)

    if doscar0.Ns == 1:
        ax[0].set_ylabel("DOS (eV⁻¹)", color=palette["black"])
        ax[0].set_xlabel("Energy (eV)", color=palette["black"], labelpad=-1)
    else:
        ax[0].set_xticklabels([])
        ax[0].set_ylabel("DOS (eV⁻¹)", color=palette["black"], y=0.0)
        ax[1].set_xlabel("Energy (eV)", color=palette["black"], labelpad=-1)

else:  # vertical
    fig = plt.figure(figsize=(1, 3.75))
    gs0 = fig.add_gridspec(1, doscar0.Ns, wspace=0.0, hspace=0.00, left=0.03, right=0.97, top=0.97, bottom=0.07)
    ax = []
    for ispin in range(doscar0.Ns):
        ax.append(fig.add_subplot(gs0[doscar0.Ns-ispin-1]))

    for ispin in range(doscar0.Ns):
        ax[ispin].axhline(linewidth=1, color=palette["gray"], zorder=0)
        ax[ispin].plot(doscar0.dos[ispin], doscar0.energy, color=palette[color_dos[ispin]], linewidth=1, zorder=3)

    if is_proj:
        for ispin in range(doscar0.Ns):
            ax[ispin].plot(pdos[ispin, :], doscar0.energy_pdos, color=palette[color_pdos[ispin]],
                           linewidth=1, zorder=3.5)
        filename_out = "dos_v_"+atom_list+"_"+orb_list+".png"

    ax[0].set_xlim([0, dosmax*1.1])
    if doscar0.Ns == 2:
        ax[1].set_xlim([dosmax*1.1, 0])
    for ispin in range(doscar0.Ns):
        ax[ispin].set_ylim([xmin, xmax])
        ax[ispin].tick_params(axis="x", bottom=False, top=False, direction="in", length=0)
        ax[ispin].tick_params(axis="y", left=True, right=True, direction="in",
                              color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0, pad=4)
        ax[ispin].set_yticklabels([])

    if doscar0.Ns == 1:
        # put "DOS" as ticklabel to align with the K point labels in bs
        ax[0].set_xticks([dosmax*0.55], ["DOS"], color=palette["black"])
    elif doscar0.Ns == 2:
        ax[0].set_xticks([0.0], ["DOS"], color=palette["black"])
        ax[1].set_xticks([], [])

for ispin in range(doscar0.Ns):
    for edge in ["bottom", "top", "left", "right"]:
        ax[ispin].spines[edge].set_color(palette["black"])
        ax[ispin].spines[edge].set_linewidth(1)
        ax[ispin].spines[edge].set_zorder(4)

fig.savefig(filename_out, dpi=1200)
