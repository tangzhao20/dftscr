import sys
import numpy as np
import re


class Procar:
    # Ns
    # Nk
    # Nb
    # Na
    # Norb
    # proj[Ns][Nk][Nb][Na][Norb]  # numpy
    # orb_name[Norb]

    def __init__(self):
        self.Ns = 1
        self.Nk = 0
        self.Nb = 0
        self.Na = 0
        self.Norb = 0
        self.proj = None
        self.orb_name = []

    def read_vasp(self, filename="PROCAR"):

        f0 = open(filename, "r")
        line = f0.readlines()
        f0.close()
        word = line[1].split()
        self.Nk = int(word[3])
        self.Nb = int(word[7])
        self.Na = int(word[11])

        if len(line) > (((self.Na+5)*self.Nb+3)*self.Nk+1)*1.5:
            # the file length should be (((Na+5)Nb+3)Nk+1)Ns+1
            self.Ns = 2

        word = line[7].split()
        self.Norb = len(word)-2

        print("Ns: "+str(self.Ns)+"  Na: "+str(self.Na)+"  Nk: "+str(self.Nk) +
              "  Nb: "+str(self.Nb)+"  Norb: "+str(self.Norb))

        self.proj = np.zeros((self.Ns, self.Nk, self.Nb, self.Na, self.Norb))
        ispin = -1
        for il in range(len(line)):
            word = line[il].split()
            if len(word) == 0 or word[0][0] == "!":
                continue
            if len(word) >= 9 and word[0] == "#" and word[2] == "k-points:":
                ispin += 1
            elif word[0] == "k-point":
                ik = int(word[1]) - 1
            elif word[0] == "band":
                ib = int(word[1]) - 1
            elif word[0].isdigit():
                ia = int(word[0]) - 1
                for iorb in range(self.Norb):
                    self.proj[ispin, ik, ib, ia, iorb] = float(word[iorb+1])

        if self.Norb == 9:
            self.orb_name = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "x2-y2"]
        elif self.Norb == 16:
            self.orb_name = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz",
                             "x2-y2", "fy3x2", "fxyz", "fyz2", "fz3", "fxz2", "fzx2", "fx3"]
        else:
            print("Norb = "+str(self.Norb)+" is not support yet")
            sys.exit()

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

    def plot(self, atom_flag, orb_flag):
        # plot_proj[Ns][Nb][Nk]  # numpy
        plot_proj = self.proj[:, :, :, atom_flag, :][:, :, :, :, orb_flag].sum(axis=(3, 4)).swapaxes(1, 2)
        return plot_proj

    def read_orb_list(self, orb_string):
        orb_list = []
        orb_flag = np.zeros(self.Norb, dtype=bool)
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
                orb_flag[self.orb_name.index(orb)] = True
            else:
                print("projector "+orb+" does not exist")

        return orb_flag
