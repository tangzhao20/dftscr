#!/usr/bin/env python3

# plot DOS

# python dos.py (v) package (E1) (E2)

# VASP input: DOSCAR
# QE input: *.dos *.xml
# QE projection input: *.dos *.xml *.pdos_atm#*(*)_wfc#*(*)

import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from classes import Doscar
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
if package in package_name["vasp"]:
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

lproj = False
if package in package_name["qeproj"]:
    lproj = True
    doscar0.readpdos_qe()
    if doscar0.Ns == 2:
        print("Spin polarized projected band structure not supported yet.\n")
        sys.exit()

    # Input: *.pdos_atm#*(*)_wfc#*(*)
    doscar0.readpdos_qe()

    atom_list = sys.argv[2]
    del sys.argv[2]
    atom_flag = doscar0.read_atom_list(atom_list)
    orb_list = sys.argv[2]
    del sys.argv[2]
    orb_flag = doscar0.read_orb_list(orb_list)

if len(sys.argv) >= 4:
    xmax = float(sys.argv[3])
    xmin = float(sys.argv[2])
elif len(sys.argv) == 3:
    xmax = float(sys.argv[2])
    xmin = -float(sys.argv[2])
else:
    xmax = 5.0
    xmin = -5.0

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

if fvertical == False:
    fig = plt.figure(figsize=(5, 3.75))
    if doscar0.Ns == 1:
        gs0 = fig.add_gridspec(1, 1, wspace=0.0, hspace=0.00, left=0.14, right=0.96, top=0.97, bottom=0.12)
        ax0 = gs0.subplots()
    else:
        gs0 = fig.add_gridspec(2, 1, wspace=0.0, hspace=0.00, left=0.14, right=0.96, top=0.97, bottom=0.12)
        (ax0, ax1) = gs0.subplots()

    ax0.axvline(linewidth=1, color=palette["gray"], zorder=0)
    ax0.plot(doscar0.energy, doscar0.dos[0], color=palette["darkblue"], linewidth=1, zorder=3)
    if package in package_name["qeproj"]:
        pdos = doscar0.plot_pdos(atom_flag, orb_flag)
        ax0.plot(doscar0.energy_pdos, pdos, color=palette["orange"], linewidth=1, zorder=3.5)
        filename_out = "dos_"+atom_list+"_"+orb_list+".png"

    ax0.set_xlim([xmin, xmax])
    ax0.set_ylim([0, dosmax*1.1])
    ax0.tick_params(axis="x", bottom=True, top=True, direction="in",
                    color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0, pad=4)

else:  # vertical
    fig = plt.figure(figsize=(1, 3.75))
    if doscar0.Ns == 1:
        gs0 = fig.add_gridspec(1, 1, wspace=0.0, hspace=0.00, left=0.03, right=0.97, top=0.97, bottom=0.07)
        ax0 = gs0.subplots()
    else:
        gs0 = fig.add_gridspec(1, 2, wspace=0.0, hspace=0.00, left=0.03, right=0.97, top=0.97, bottom=0.07)
        (ax1, ax0) = gs0.subplots()

    ax0.axhline(linewidth=1, color=palette["gray"], zorder=0)
    ax0.plot(doscar0.dos[0], doscar0.energy, color=palette["darkblue"], linewidth=1, zorder=3)
    if package in package_name["qeproj"]:
        pdos = doscar0.plot_pdos(atom_flag, orb_flag)
        ax0.plot(pdos, doscar0.energy_pdos, color=palette["orange"], linewidth=1, zorder=3.5)
        filename_out = "dos_v_"+atom_list+"_"+orb_list+".png"

    ax0.set_xlim([0, dosmax*1.1])
    ax0.set_ylim([xmin, xmax])
    ax0.tick_params(axis="x", bottom=False, top=False, direction="in", length=0)

ax0.tick_params(axis="y", left=True, right=True, direction="in",
                color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0, pad=4)
for edge in ["bottom", "top", "left", "right"]:
    ax0.spines[edge].set_color(palette["black"])
    ax0.spines[edge].set_linewidth(1)
    ax0.spines[edge].set_zorder(4)

if fvertical == False:
    if doscar0.Ns == 1:
        ax0.set_xlabel("Energy (eV)", color=palette["black"], labelpad=-1)
        ax0.set_ylabel("DOS (eV⁻¹)", color=palette["black"])

    elif doscar0.Ns == 2:
        ax0.set_xticklabels([])
        ax1.axvline(linewidth=1, color=palette["gray"], zorder=0)
        ax1.plot(doscar0.energy, doscar0.dos[1], color=palette["orange"], linewidth=1, zorder=3)
        ax1.set_xlim([xmin, xmax])
        ax1.set_ylim([dosmax*1.1, 0.0])
        ax1.set_xlabel("Energy (eV)", color=palette["black"], labelpad=-1)
        ax0.set_ylabel("DOS (eV⁻¹)", color=palette["black"], y=0.0)
        ax1.tick_params(axis="x", bottom=True, top=True, direction="in",
                        color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0, pad=4)
else:  # vertical
    ax0.set_yticklabels([])
    if doscar0.Ns == 1:
        # put "DOS" as ticklabel to align with the K point labels in bs
        ax0.set_xticks([dosmax*0.55], ["DOS"], color=palette["black"])

    elif doscar0.Ns == 2:
        ax0.set_xticks([0.0], ["DOS"], color=palette["black"])
        ax1.set_xticks([], [])
        ax1.set_yticklabels([])
        ax1.axhline(linewidth=1, color=palette["gray"], zorder=0)
        ax1.plot(doscar0.dos[1], doscar0.energy, color=palette["orange"], linewidth=1, zorder=3)
        ax1.set_xlim([dosmax*1.1, 0.0])
        ax1.set_ylim([xmin, xmax])
        ax1.tick_params(axis="x", bottom=False, top=False, direction="in", length=0)

if doscar0.Ns == 2:
    ax1.tick_params(axis="y", left=True, right=True, direction="in",
                    color=palette["gray"], labelcolor=palette["black"], width=1, zorder=0, pad=4)
    for edge in ["bottom", "top", "left", "right"]:
        ax1.spines[edge].set_color(palette["black"])
        ax1.spines[edge].set_linewidth(1)
        ax1.spines[edge].set_zorder(4)

fig.savefig(filename_out, dpi=1200)
