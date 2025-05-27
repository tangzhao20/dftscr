import sys
import numpy as np
import re


class Procar:
    # Nk
    # Nb
    # Na
    # Norb
    # proj[Nk][Nb][Na][Norb]
    # orb_name[Norb]

    def __init__(self):
        self.Nk = 0
        self.Nb = 0
        self.Na = 0
        self.Norb = 0
        self.proj = []
        self.orb_name = []

    def read_vasp(self, filename="PROCAR"):

        f0 = open(filename, "r")
        line = f0.readlines()
        f0.close()
        word = line[1].split()
        self.Nk = int(word[3])
        self.Nb = int(word[7])
        self.Na = int(word[11])

        word = line[7].split()
        self.Norb = len(word)-2

        print("Na "+str(self.Na)+" Nk "+str(self.Nk)+" Nb "+str(self.Nb)+" Norb "+str(self.Norb))

        self.proj = []
        ia = -1
        for il in range(len(line)):
            word = line[il].split()
            if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                continue
            if word[0] == "k-point":
                ik = int(word[1])
                proj0 = []
                # word=line[].split()
            elif word[0] == "band":
                proj1 = []
                ib = int(word[1])
            elif word[0] == "ion":
                ia = 0
            elif ia >= 0 and ia < self.Na:
                if int(word[0]) != ia+1:
                    print("atom number mismatch on line "+str(il+1))
                    sys.exit()
                proj2 = []
                for lm in range(self.Norb):
                    proj2.append(float(word[lm+1]))
                proj1.append(proj2)
                if ia < self.Na-1:
                    ia = ia+1
                else:
                    ia = -1
                    proj0.append(proj1)
                    if ib == self.Nb-1:
                        self.proj.append(proj0)

        if self.Norb == 9:
            self.orb_name = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "x2-y2"]
        elif self.Norb == 16:
            self.orb_name = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz",
                             "x2-y2", "fy3x2", "fxyz", "fyz2", "fz3", "fxz2", "fzx2", "fx3"]
        else:
            print("Norb = "+str(self.Norb)+" is not support yet")
            sys.exit()

        del proj0
        del proj1
        del proj2
        del line
        del word

    def read_qe(self, filename="projwfc.out"):
        f0 = open(filename, "r")
        line = f0.readlines()
        f0.close()

        lmmap = []
        atommap = []
        fread = False
        ffirstk = True
        for il in range(len(line)):
            word = line[il].split()
            if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                continue
            if word[0] == "state":
                ia = int(word[4])-1
                if len(atommap) == 0 or ia != atommap[-1]:
                    iorb = 0
                else:
                    iorb += 1
                lmmap.append(iorb)
                atommap.append(ia)
                self.Norb = max(iorb+1, self.Norb)
                self.Na = max(ia+1, self.Na)
            elif word[0] == "nkstot":
                self.Nk = int(word[2])
            elif word[0] == "nbnd":
                self.Nb = int(word[2])
            elif word[0] == "k":
                if ffirstk:
                    ffirstk = False
                    for ik in range(self.Nk):
                        proj0 = []
                        for ib in range(self.Nb):
                            proj1 = []
                            for ia in range(self.Na):
                                proj1.append([0.0]*self.Norb)
                            proj0.append(proj1)
                            del proj1
                        self.proj.append(proj0)
                        del proj0
                    ik = -1
                ik += 1
                ib = -1
            elif word[0] == "psi":
                ib += 1
                fread = True
            elif word[0] == "|psi|^2":
                fread = False

            if fread:
                number = re.findall(r"[-+]?(?:\d*\.*\d+)", line[il])  # find all numbers
                for ii in range(0, len(number), 2):
                    # proj[k][b][a][orb]
                    self.proj[ik][ib][atommap[int(number[ii+1])-1]][lmmap[int(number[ii+1])-1]] = float(number[ii])
        if self.Norb == 4:
            self.orb_name = ["s", "pz", "px", "py"]
        elif self.Norb == 9:
            self.orb_name = ["s", "pz", "px", "py", "dz2", "dxz", "dyz", "x2-y2", "dxy"]
        elif self.Norb == 16:
            self.orb_name = ["s", "pz", "px", "py", "dz2", "dxz", "dyz", "x2-y2",
                             "dxy", "fz3", "fxz2", "fyz2", "fzx2", "fxyz", "fx3", "fy3x2"]
        else:
            print("Norb = "+str(self.Norb)+" is not support yet")
            sys.exit()

    def plot(self, atom_flag, orb_flag, fac):
        plot_proj = []
        for ib in range(self.Nb):
            plot_proj0 = []
            for ik in range(self.Nk):
                plot_proj1 = 0.0
                for ia in range(self.Na):
                    for iorb in range(self.Norb):
                        if atom_flag[ia] and orb_flag[iorb] == 1:
                            plot_proj1 += self.proj[ik][ib][ia][iorb]
                plot_proj0.append(plot_proj1*fac)
            plot_proj.append(plot_proj0)
        return plot_proj

    def read_orb_list(self, orb_string):
        orb_list = []
        orb_flag = [0]*self.Norb
        for orb in orb_string.split("+"):
            if orb in self.orb_name:
                orb_list.append(orb)
            elif orb == "p":
                orb_list += ["px", "py", "pz"]
            elif orb == "d":
                orb_list += ["dxy", "dyz", "dz2", "dxz", "x2-y2"]
            elif orb == "f":
                orb_list += ["fy3x2", "fxyz", "fyz2", "fz3", "fxz2", "fzx2", "fx3"]
            elif orb == "dx2-y2":
                orb_list.append("x2-y2")
            elif orb == "all":
                orb_list += self.orb_name
            else:
                print("projector "+orb+" does not exist")

        for orb in orb_list:
            if orb in self.orb_name:
                orb_flag[self.orb_name.index(orb)] = 1
            else:
                print("projector "+orb+" does not exist")

        return orb_flag

    def read_atom_list(self, atom_string, poscar1):
        atom_flag = []
        if atom_string == "all":
            for ia in range(self.Na):
                atom_flag.append(1)
        else:
            for ia in range(self.Na):
                atom_flag.append(0)
            for atom_section in atom_string.split(","):
                is_valid_atom = False
                ia_start = 0
                ia_end = 0
                for itype in range(len(poscar1.atomtype)):
                    ia_end = ia_end + poscar1.Naint[itype]
                    if atom_section == poscar1.atomtype[itype]:
                        for ia in range(ia_start, ia_end):
                            atom_flag[ia] = 1
                        is_valid_atom = True
                    ia_start = ia_end
                if is_valid_atom == False:
                    atom_list = atom_section.split("-")
                    if len(atom_list) == 1:
                        atom_flag[int(atom_list[0])-1] = 1
                        is_valid_atom = True
                    else:
                        for ia in range(int(atom_list[0])-1, int(atom_list[1])):
                            atom_flag[ia] = 1
                        is_valid_atom = True
        return atom_flag
