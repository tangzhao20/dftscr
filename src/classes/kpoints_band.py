import sys
import os


class KpointsBand:
    """
    A class to handle k-points along high symmetry lines for band structure calculations.

    Attributes:
        kp_dict (dict): Mapping k-point labels to their fractional coordinates.
        x_labels[ip][ik]: Containing paths of high symmetry k-point labels.
        Nk_line[ip][ik]: A list of lists containing the number of k-points along each segment in the paths.

    # Use kp_dict[x_labels[ip][ik]][ix] for the ik-th points in the ip-th path. ix indicates kx, ky, or kz
    """

    def __init__(self):
        self.kp_dict = {}
        self.x_labels = []
        self.Nk_line = []

    def __str__(self):
        str_out = "KPOINTS:\n"
        str_out += " Nk_line = " + str(self.Nk_line) + "\n"
        for x_label in self.x_labels:
            str_out += " "
            for x_label0 in x_label:
                str_out += " " + x_label0
            str_out += "\n"
        return str_out

#########################################################################

    def read_vasp(self, filename="KPOINTS"):

        f0 = open(filename, "r")
        line = f0.readlines()
        f0.close()
        if line[2][0] not in {"L", "l"}:
            print("line-mode required.")
            sys.exit()
        if line[3][0] in {"C", "c", "K", "k"}:
            print("fractional coordinates of k-points required.")
            sys.exit()

        Nk_line = int(line[1].split()[0])  # kpoints on a line

        greek = self.loadgreek()

        is_first = True
        for il in range(4, len(line)):
            word = line[il].split()
            if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                continue
            word = line[il].replace("!", " ").replace("#", " ").split()
            label = word[3]
            if label in greek:
                label = greek[label]
            if label not in self.kp_dict:
                self.kp_dict[label] = [float(word[0]), float(word[1]), float(word[2])]
            if is_first:
                if len(self.x_labels) == 0 or label != self.x_labels[-1][-1]:
                    self.x_labels.append([label])
                    self.Nk_line.append([Nk_line])
            else:
                self.x_labels[-1].append(label)
                self.Nk_line[-1].append(Nk_line)
            is_first = not is_first

    def read_kpathin(self, filename="kpath.in", Nk_line=40):
        # self.Nk_line=Nk_line

        f0 = open(filename, "r")
        line = f0.readlines()
        f0.close()

        greek = self.loadgreek()

        for l in line:
            word = l.split()
            if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                continue
            # if len(word)==4 :
            try:
                float(word[0])
                float(word[1])
                float(word[2])
            except:
                self.x_labels.append(word)
                self.Nk_line.append([Nk_line]*(len(word)-1))
                for ik in range(len(self.x_labels[-1])):
                    if self.x_labels[-1][ik] in greek:
                        self.x_labels[-1][ik] = greek[self.x_labels[-1][ik]]
            else:
                label = word[3]
                if label in greek:
                    label = greek[label]
                self.kp_dict[label] = [float(word[0]), float(word[1]), float(word[2])]

#########################################################################

    def write_vasp(self, filename="KPOINTS.new"):
        f1 = open(filename, "w")
        f1.write("k points along high symmetry lines\n")
        f1.write(f"{self.Nk_line[0][0]:d}\n")
        f1.write("line mode\n")
        f1.write("fractional\n")
        for p in self.x_labels:
            for ih in range(len(p)-1):
                for ix in range(3):
                    f1.write("{:16.12f}".format(self.kp_dict[p[ih]][ix])+" ")
                f1.write(p[ih]+"\n")
                for ix in range(3):
                    f1.write("{:16.12f}".format(self.kp_dict[p[ih+1]][ix])+" ")
                f1.write(p[ih+1]+"\n\n")
        f1.close()

    def write_qe(self, filename="kpath.out", Nk=0, rlc=None):
        if Nk > 0:
            x_ticks = self.x_ticks_out(rlc)

            x_total = 0.0
            for ip in range(len(x_ticks)):
                x_total += x_ticks[ip][-1]
            for ip in range(len(x_ticks)):
                for ih in range(1, len(x_ticks[ip])):
                    self.Nk_line[ip][ih-1] = int((x_ticks[ip][ih]-x_ticks[ip][ih-1]) / x_total * Nk)

            Nktotal = 0
            for ip in range(len(x_ticks)):
                Nktotal += sum(self.Nk_line[ip])

        kplot = []
        for ip in range(len(self.x_labels)):
            kplot.append(self.kp_dict[self.x_labels[ip][0]])
            for ih in range(len(self.x_labels[ip])-1):
                for ik in range(self.Nk_line[ip][ih]):
                    k = [0.0, 0.0, 0.0]
                    for ix in range(3):
                        k[ix] = (self.kp_dict[self.x_labels[ip][ih]][ix] *
                                 float(self.Nk_line[ip][ih]-ik-1) / float(self.Nk_line[ip][ih]) +
                                 self.kp_dict[self.x_labels[ip][ih+1]][ix] *
                                 float(ik+1) / float(self.Nk_line[ip][ih]))
                    kplot.append(k)

        f1 = open(filename, "w")
        f1.write("K_POINTS crystal_b\n")
        f1.write(str(len(kplot))+"\n")
        for ik in range(len(kplot)):
            for ix in range(3):
                f1.write("{:16.12f}".format(kplot[ik][ix])+" ")
            f1.write("1.0\n")
        f1.close()

    def write_kpathin(self, filename="kpath.in.new"):
        f1 = open(filename, "w")
        for h in self.kp_dict:
            for ix in range(3):
                f1.write(str(self.kp_dict[h][ix])+" ")
            f1.write(h+"\n")
        f1.write("\n")
        for p in self.x_labels:
            for h in p:
                f1.write(h+" ")
            f1.write("\n")
        f1.close()

    def write_wannier90(self, filename="wannier90_kpath.dat"):
        f1 = open(filename, "w")
        f1.write("bands_plot = true\n\n")
        f1.write("begin kpoint_path\n")
        for p in self.x_labels:
            for ih in range(len(p)-1):
                f1.write(p[ih]+" ")
                for ix in range(3):
                    f1.write(f"{self.kp_dict[p[ih]][ix]:12.8f}")
                f1.write(" "+p[ih+1]+" ")
                for ix in range(3):
                    f1.write(f"{self.kp_dict[p[ih+1]][ix]:12.8f}")
                f1.write("\n")
        f1.write("end kpoint_path\n\n")
        f1.close()

    def write_parsec(self, filename="parsec_kpath.dat"):
        f1 = open(filename, "w")
        f1.write("begin bandstruc\n")
        ih2 = 1
        for p in self.x_labels:
            for ih in range(len(p)-1):
                f1.write(f"{ih2:4d}")
                for ix in range(3):
                    f1.write(f"{self.kp_dict[p[ih]][ix]:12.8f}")
                for ix in range(3):
                    f1.write(f"{self.kp_dict[p[ih+1]][ix]:12.8f}")
                f1.write(" "+p[ih]+"-"+p[ih+1]+"\n")
                ih2 += 1
        f1.write("end bandstruc\n")
        f1.write("bandstruc_points 40\n\n")
        f1.close()

#########################################################################

    def x_labels_out(self):
        x_labels = [self.x_labels[0][:]]
        for ip in range(1, len(self.x_labels)):
            x_labels[-1][-1] = x_labels[-1][-1]+"|"+self.x_labels[ip][0]
            x_labels.append([""])
            x_labels[-1] = x_labels[-1]+self.x_labels[ip][1:]
        return x_labels

    def x_ticks_out(self, rlc):
        x_ticks = []
        for p in self.x_labels:
            x_ticks.append([0.0])
            for ik in range(len(p)):
                if ik > 0:
                    kpcold = [kpc[0], kpc[1], kpc[2]]
                    del kpc
                kpc = []
                for i in range(3):
                    kpc0 = 0.0
                    for j in range(3):
                        kpc0 = kpc0+self.kp_dict[p[ik]][j]*rlc[i, j]
                    kpc.append(kpc0)
                if ik > 0:
                    dkpc = ((kpc[0]-kpcold[0])**2+(kpc[1]-kpcold[1])**2+(kpc[2]-kpcold[2])**2)**0.5
                    x_ticks[-1].append(x_ticks[-1][-1]+dkpc)
                    del kpcold
        return (x_ticks)

    def loadgreek(self):
        this_dir, this_filename = os.path.split(__file__)
        DATA_PATH = os.path.join(this_dir, "..", "..", "data", "greek.dat")
        f0 = open(DATA_PATH, "r")
        line = f0.readlines()
        f0.close()
        greek = {}
        for l in line:
            word = l.split()
            if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                continue
            if len(word) != 2:
                print("Error: check the structure of dftscr/data/greek.dat")
                sys.exit()
            greek[word[0]] = word[1]
        return greek

    def findlabel(self, k, dim=0):
        # k[3] is a k point in direct coordinate,
        # this function find the label of k
        # dim is the maximum dimension of searching. (point, line, plane, ...)
        tol = 1e-4
        if dim >= 0:
            for h in self.kp_dict:
                kh = self.kp_dict[h]
                if abs(kh[0]-k[0]) <= tol and abs(kh[1]-k[1]) <= tol and abs(kh[2]-k[2]) <= tol:
                    return h
        if dim >= 1:
            for h1 in self.kp_dict:
                kh1 = self.kp_dict[h1]
                for h2 in self.kp_dict:
                    if h1 == h2:
                        continue
                    kh2 = self.kp_dict[h2]
                    x1 = kh1[0]-k[0]
                    x2 = kh2[0]-k[0]
                    y1 = kh1[1]-k[1]
                    y2 = kh2[1]-k[1]
                    z1 = kh1[2]-k[2]
                    z2 = kh2[2]-k[2]
                    if abs(x1*y2-x2*y1) <= tol and abs(y1*z2-y2*z1) <= tol and abs(x1*z2-x2*z1) <= tol:
                        return "("+h1+","+h2+")"
        return "elsewhere"
