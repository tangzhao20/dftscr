import sys
import os


class KpointsBand:
    def __init__(self):

        self.kpdict = {}  # a dictionary which maps label to b-coordinate
        self.xlabels = []  # a path of high symmetry k-points. xlabels[Np][Nh[ip]]
        # use kpdict[xlabels[ip][ik]][ix] for the ik-th points in the ip-th path. ix indicates kx, ky, or kz
        self.Nk_line = []

    def __str__(self):
        str_out = "KPOINTS:\n"
        str_out += " Nk_line = " + str(self.Nk_line) + "\n"
        for xlabel in self.xlabels:
            str_out += " "
            for xlabel0 in xlabel:
                str_out += " " + xlabel0
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

        ip = 0
        ffirst = True
        for il in range(4, len(line)):
            word = line[il].split()
            if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                continue
            word = line[il].replace("!", " ").replace("#", " ").split()
            label = word[3]
            if label in greek:
                label = greek[label]
            if label not in self.kpdict:
                self.kpdict[label] = [float(word[0]), float(word[1]), float(word[2])]
            if ffirst == True:
                if len(self.xlabels) == 0 or label != self.xlabels[-1][-1]:
                    self.xlabels.append([label])
                    self.Nk_line.append([Nk_line])
            else:
                self.xlabels[-1].append(label)
                self.Nk_line[-1].append(Nk_line)
            ffirst = not ffirst

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
                self.xlabels.append(word)
                self.Nk_line.append([Nk_line]*(len(word)-1))
                for ik in range(len(self.xlabels[-1])):
                    if self.xlabels[-1][ik] in greek:
                        self.xlabels[-1][ik] = greek[self.xlabels[-1][ik]]
            else:
                label = word[3]
                if label in greek:
                    label = greek[label]
                self.kpdict[label] = [float(word[0]), float(word[1]), float(word[2])]

#########################################################################

    def write_vasp(self, filename="KPOINTS.new"):
        f1 = open(filename, "w")
        f1.write("k points along high symmetry lines\n")
        f1.write("40\n")
        f1.write("line mode\n")
        f1.write("fractional\n")
        for p in self.xlabels:
            for ih in range(len(p)-1):
                for ix in range(3):
                    f1.write("{:16.12f}".format(self.kpdict[p[ih]][ix])+" ")
                f1.write(p[ih]+"\n")
                for ix in range(3):
                    f1.write("{:16.12f}".format(self.kpdict[p[ih+1]][ix])+" ")
                f1.write(p[ih+1]+"\n\n")
        f1.close()

    def write_qe(self, filename="kpath.out", Nk=0, rlc=[]):
        if Nk > 0:
            xticks = self.xticks_out(rlc=rlc)

            Ltotal = 0.0
            for ip in range(len(xticks)):
                Ltotal += xticks[ip][-1]
            for ip in range(len(xticks)):
                for ih in range(1, len(xticks[ip])):
                    self.Nk_line[ip][ih-1] = int((xticks[ip][ih]-xticks[ip][ih-1])/Ltotal*Nk)

            Nktotal = 0
            for ip in range(len(xticks)):
                Nktotal += sum(self.Nk_line[ip])

        kplot = []
        for ip in range(len(self.xlabels)):
            kplot.append(self.kpdict[self.xlabels[ip][0]])
            for ih in range(len(self.xlabels[ip])-1):
                for ik in range(self.Nk_line[ip][ih]):
                    k = [0.0, 0.0, 0.0]
                    for ix in range(3):
                        k[ix] = self.kpdict[self.xlabels[ip][ih]][ix]*float(self.Nk_line[ip][ih]-ik-1)/float(
                            self.Nk_line[ip][ih])+self.kpdict[self.xlabels[ip][ih+1]][ix]*float(ik+1)/float(self.Nk_line[ip][ih])
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
        for h in self.kpdict:
            for ix in range(3):
                f1.write(str(self.kpdict[h][ix])+" ")
            f1.write(h+"\n")
        f1.write("\n")
        for p in self.xlabels:
            for h in p:
                f1.write(h+" ")
            f1.write("\n")
        f1.close()

    def write_wannier90(self, filename="wannier90_kpath.dat"):
        f1 = open(filename, "w")
        f1.write("bands_plot = true\n\n")
        f1.write("begin kpoint_path\n")
        for p in self.xlabels:
            for ih in range(len(p)-1):
                f1.write(p[ih]+" ")
                for ix in range(3):
                    f1.write(f"{self.kpdict[p[ih]][ix]:12.8f}")
                f1.write(" "+p[ih+1]+" ")
                for ix in range(3):
                    f1.write(f"{self.kpdict[p[ih+1]][ix]:12.8f}")
                f1.write("\n")
        f1.write("end kpoint_path\n\n")
        f1.close()

    def write_parsec(self, filename="parsec_kpath.dat"):
        f1 = open(filename, "w")
        f1.write("begin bandstruc\n")
        ih2 = 1
        for p in self.xlabels:
            for ih in range(len(p)-1):
                f1.write(f"{ih2:4d}")
                for ix in range(3):
                    f1.write(f"{self.kpdict[p[ih]][ix]:12.8f}")
                for ix in range(3):
                    f1.write(f"{self.kpdict[p[ih+1]][ix]:12.8f}")
                f1.write(" "+p[ih]+"-"+p[ih+1]+"\n")
                ih2 += 1
        f1.write("end bandstruc\n")
        f1.write("bandstruc_points 40\n\n")
        f1.close()

#########################################################################

    def xlabels_out(self):
        xlabels = [self.xlabels[0][:]]
        for ip in range(1, len(self.xlabels)):
            xlabels[-1][-1] = xlabels[-1][-1]+"|"+self.xlabels[ip][0]
            xlabels.append([""])
            xlabels[-1] = xlabels[-1]+self.xlabels[ip][1:]
        return xlabels

    def xticks_out(self, rlc):
        xticks = []
        for p in self.xlabels:
            xticks.append([0.0])
            for ik in range(len(p)):
                if ik > 0:
                    kpcold = [kpc[0], kpc[1], kpc[2]]
                    del kpc
                kpc = []
                for i in range(3):
                    kpc0 = 0.0
                    for j in range(3):
                        kpc0 = kpc0+self.kpdict[p[ik]][j]*rlc[i][j]
                    kpc.append(kpc0)
                if ik > 0:
                    dkpc = ((kpc[0]-kpcold[0])**2+(kpc[1]-kpcold[1])**2+(kpc[2]-kpcold[2])**2)**0.5
                    xticks[-1].append(xticks[-1][-1]+dkpc)
                    del kpcold
        return (xticks)

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
            for h in self.kpdict:
                kh = self.kpdict[h]
                if abs(kh[0]-k[0]) <= tol and abs(kh[1]-k[1]) <= tol and abs(kh[2]-k[2]) <= tol:
                    return h
        if dim >= 1:
            for h1 in self.kpdict:
                kh1 = self.kpdict[h1]
                for h2 in self.kpdict:
                    if h1 == h2:
                        continue
                    kh2 = self.kpdict[h2]
                    x1 = kh1[0]-k[0]
                    x2 = kh2[0]-k[0]
                    y1 = kh1[1]-k[1]
                    y2 = kh2[1]-k[1]
                    z1 = kh1[2]-k[2]
                    z2 = kh2[2]-k[2]
                    if abs(x1*y2-x2*y1) <= tol and abs(y1*z2-y2*z1) <= tol and abs(x1*z2-x2*z1) <= tol:
                        return "("+h1+","+h2+")"
        return "elsewhere"
