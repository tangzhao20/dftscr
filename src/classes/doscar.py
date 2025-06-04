import sys
import os
import re
import numpy as np
import xml.etree.ElementTree as ET
from load_data import load_constant


class Doscar:
    # ef : float
    # Ns : integer
    # Nedos : integer
    # energy[Nedos]
    # dos[Ns][Nedos]

    def __init__(self):
        self.ef = 0.0
        self.Ns = 0
        self.Nedos = 0
        self.energy = []
        self.dos = [[]]
        self.has_pdos = False

    def __str__(self):
        str_out = "DOSCAR:\n"
        str_out += " Ns = " + str(self.Ns) + "\n"
        str_out += " Nedos = " + str(self.Nedos) + "\n"
        return str_out

#######################################################################

    def read_vasp(self, filename="DOSCAR"):
        f0 = open(filename, "r")
        line = f0.readlines()
        f0.close()
        if len(line) < 7:
            print("Error: incomplete DOSCAR")
            sys.exit()
        self.Nedos = int(line[5].split()[2])
        self.ef = float(line[5].split()[3])
        self.Ns = (len(line[6].split()) - 1)//2
        self.energy = [0.0]*self.Nedos
        if self.Ns == 1:
            self.dos = [[0.0]*self.Nedos]
        else:
            self.dos = [[0.0]*self.Nedos, [0.0]*self.Nedos]

        for il in range(self.Nedos):
            word = line[il + 6].split()
            self.energy[il] = float(word[0])
            self.dos[0][il] = float(word[1])
            if self.Ns == 2:
                self.dos[1][il] = float(word[2])

    def read_xml(self, filename=""):
        # read dos from .xml file
        Ha = load_constant("rydberg")*2.0

        if filename == "":
            # find a .xml file
            files = os.listdir()
            for f in files:
                if f.endswith('.xml'):
                    filename = f
                    break
        tree = ET.parse(filename)
        # TODO: test the metal case
        if tree.getroot().find("output").find("band_structure").find("highestOccupiedLevel") is None:
            self.ef = float(tree.getroot().find("output").find("band_structure").find("fermi_energy").text)*Ha
        else:
            self.ef = float(tree.getroot().find("output").find("band_structure").find("highestOccupiedLevel").text)*Ha

    def read_qe(self, filename="pwscf.dos"):
        f0 = open(filename, "r")
        line = f0.readlines()
        f0.close()

        if len(line) < 2:
            print("Error: incomplete " + filename)
            sys.exit()
        self.ef = float(line[0].split()[-2])
        self.Ns = (len(line[1].split()) - 1)//2
        if self.Ns == 2:
            self.dos.append([])
        for l in line:
            word = l.split()
            if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                continue
            self.Nedos += 1
            self.energy.append(float(word[0]))
            self.dos[0].append(float(word[1]))
            if self.Ns == 2:
                self.dos[1].append(float(word[2]))

    def readpdos_qe(self):
        # pdos[Ns][Nepdos][Na][Norb]
        # x_pdos[Nepdos]
        # atomtype[Ntype]
        # Naint[Ntype]
        # orb_name[Norb]

        if self.Ns == 2:
            print("QE pdos has not been supported for spin polarized system yet.")
            sys.exit()

        files = os.listdir()
        file_list = []
        for f in files:
            fileinfo = re.search(r"(.+?)\.pdos_atm#(\d+)\((.+?)\)_wfc#(\d+)\((.+?)\)", f)
            if fileinfo:
                filename = fileinfo.group(0)
                system = fileinfo.group(1)
                iatom = fileinfo.group(2)
                atomtype = fileinfo.group(3)
                iorb = fileinfo.group(4)
                orb_info = fileinfo.group(5)
                file_list.append([filename, int(iatom), atomtype, int(iorb)])
        if len(file_list) == 0:
            print("Error: No valid pdos file found.")
            sys.exit()
        file_list = sorted(file_list, key=lambda x: (x[1], x[3]))
        self.Na = max(file_list, key=lambda x: x[1])[1]
        self.Norb = max(file_list, key=lambda x: x[3])[3]**2

        # Read energies and Nepdos from the first file
        f = file_list[0]
        f0 = open(f[0], "r")
        line = f0.readlines()
        f0.close()
        self.energy_pdos = []
        for l in line:
            word = l.split()
            if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                continue
            self.energy_pdos.append(float(word[0]))
        self.Nepdos = len(self.energy_pdos)

        self.pdos = []
        self.pdos = np.arange(self.Ns*self.Nepdos*self.Na*self.Norb)
        self.pdos = np.reshape(self.pdos, (self.Ns, self.Nepdos, self.Na, self.Norb)).tolist()

        for f in file_list:
            # f=[filename, iatom, atomtype, iorb]
            f0 = open(f[0], "r")
            line = f0.readlines()
            f0.close()
            ie = 0
            for l in line:
                word = l.split()
                if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                    continue
                if f[3] == 1:
                    self.pdos[0][ie][f[1] - 1][0] = float(word[2])
                elif f[3] == 2:
                    self.pdos[0][ie][f[1] - 1][1:4] = [float(word[2]), float(word[3]), float(word[4])]
                ie += 1

        # build atom list atomtype[Ntype] and Naint[Ntype]
        self.atomtype = []
        self.Naint = []
        for f in file_list:
            if f[3] == 1:
                if len(self.atomtype) == 0 or self.atomtype[-1] != f[2]:
                    self.atomtype.append(f[2])
                    self.Naint.append(1)
                else:
                    self.Naint[-1] += 1

        # Only s and p orbitals are supported now
        self.orb_name = ["s", "pz", "px", "py"]
        self.has_pdos = True

#######################################################################

    def energyshift(self, ezero):
        for ie in range(self.Nedos):
            self.energy[ie] = self.energy[ie] - ezero
        if self.has_pdos:
            for ie in range(self.Nepdos):
                self.energy_pdos[ie] = self.energy_pdos[ie] - ezero

    def plot(self, atom_flag, orb_flag):
        # plot_proj[Ns][Ne]  # numpy
        plot_pdos = self.pdos[:, :, atom_flag, :][:, :, :, orb_flag].sum(axis=(2, 3))
        return plot_pdos

        plot_pdos = []
        for ie in range(self.Nepdos):
            pdos_out0 = 0.0
            for ia in range(self.Na):
                for im in range(self.Norb):
                    if atom_flag[ia] and orb_flag[im] == 1:
                        pdos_out0 += self.pdos[0][ie][ia][im]
            plot_pdos.append(pdos_out0)
        return plot_pdos

    def read_orb_list(self, orb_list):
        orb_list0 = orb_list.split("+")
        orb_list1 = []
        orb_flag = [0]*self.Norb
        for j in orb_list0:
            if j in self.orb_name:
                orb_list1.append(j)
            elif j == "p":
                orb_list1 += ["px", "py", "pz"]
            elif j == "d":
                orb_list1 += ["dxy", "dyz", "dz2", "dxz", "x2-y2"]
            elif j == "f":
                orb_list1 += ["fy3x2", "fxyz", "fyz2", "fz3", "fxz2", "fzx2", "fx3"]
            elif j == "dx2-y2":
                orb_list1.append("x2-y2")
            elif j == "all":
                orb_list1 += self.orb_name
            else:
                print("projector " + j + " does not exist")

        for j in orb_list1:
            if j in self.orb_name:
                orb_flag[self.orb_name.index(j)] = 1
            else:
                print("projector " + j + " does not exist")

        return orb_flag

    def read_atom_list(self, atom_list):
        atom_flag = []
        if atom_list == "all":
            for a in range(self.Na):
                atom_flag.append(1)
        else:
            for a in range(self.Na):
                atom_flag.append(0)
            atom_list0 = atom_list.split(",")
            for i in range(len(atom_list0)):
                Fatom_valid = False
                Iatom0 = 0
                Iatom1 = 0
                for j in range(len(self.atomtype)):
                    Iatom1 = Iatom1 + self.Naint[j]
                    if atom_list0[i] == self.atomtype[j]:
                        for a in range(Iatom0, Iatom1):
                            atom_flag[a] = 1
                        Fatom_valid = True
                    Iatom0 = Iatom1
                if not Fatom_valid:
                    atom_list1 = atom_list0[i].split("-")
                    if len(atom_list1) == 1:
                        atom_flag[int(atom_list1[0])-1] = 1
                        Fatom_valid = True
                    else:
                        for j in range(int(atom_list1[0])-1, int(atom_list1[1])):
                            atom_flag[j] = 1
                        Fatom_valid = True
        return atom_flag
