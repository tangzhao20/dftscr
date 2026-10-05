import sys
import os
import numpy as np
import xml.etree.ElementTree as ET
from load_data import load_constant


class Eigenval:
    """
    Parses and analyzes band structure eigenvalue data.

    Attributes:
        Ns (int): The number of spins.
        Nk (int): The number of k-points.
        Nb (int): The number of bands.
        kp[Nk, 3] (numpy.ndarray): The k-points coordinates.
        weight[Nk] (numpy.ndarray): The weight of each k-point.
        eig[Ns, Nk, Nb] (numpy.ndarray): The eigenvalues (band energies) for each state.
        occ[Ns, Nk, Nb] (numpy.ndarray): The occupation numbers for each state.
        is_semic (bool): Whether the system is a semiconductor with a gap.
        Nvb[Ns] (list of int): Number of valence bands for each spin.
        vbm (float): Valence band maximum energy.
        cbm (float): Conduction band minimum energy.
        edg (float): Direct band gap.
        eindg (float): Indirect band gap.
    """

    def __init__(self):
        self.Ns = 1
        self.Nk = 0
        self.Nb = 0
        self.weight = None
        self.kp = None
        self.eig = None
        self.occ = None
        self.is_semic = False
        self.Nvb = None

    def __str__(self):
        str_out = "EIGENVAL:\n"
        str_out += " Nk = " + str(self.Nk) + "\n"
        str_out += " Nb = " + str(self.Nb) + "\n"
        str_out += " Ns = " + str(self.Ns) + "\n"
        if self.is_semic:
            str_out += " is_semic = True\n"
            str_out += f" vbm = {self.vbm:.2f} eV\n"
            str_out += f" cbm = {self.cbm:.2f} eV\n"
            str_out += f" Egdir = {self.edg:.2f} eV\n"
            if self.eindg < self.edg:
                str_out += f" Egind = {self.eindg:.2f} eV\n"
        else:
            str_out += " is_semic = False\n"
        return str_out

#########################################################################

    def read_vasp(self, filename="EIGENVAL", is_hse=False):

        with open(filename, "r") as f0:
            lines = f0.readlines()

        self.Ns = int(lines[0].split()[3])
        self.Nk = int(lines[5].split()[1])
        self.Nb = int(lines[5].split()[2])

        self.kp = np.zeros((self.Nk, 3))
        self.weight = np.zeros(self.Nk)
        self.eig = np.zeros((self.Ns, self.Nk, self.Nb))
        self.occ = np.zeros((self.Ns, self.Nk, self.Nb))

        for ik in range(self.Nk):
            offset = 7 + ik * (self.Nb + 2)
            word = lines[offset].split()
            self.kp[ik] = [float(word[0]), float(word[1]), float(word[2])]
            self.weight[ik] = float(word[3])
            for ib in range(self.Nb):
                word = lines[offset + 1 + ib].split()
                if self.Ns == 1:
                    self.eig[0, ik, ib] = float(word[1])
                    self.occ[0, ik, ib] = float(word[2])
                elif self.Ns == 2:
                    self.eig[0, ik, ib] = float(word[1])
                    self.eig[1, ik, ib] = float(word[2])
                    self.occ[0, ik, ib] = float(word[3])
                    self.occ[1, ik, ib] = float(word[4])

        if is_hse:
            mask = np.abs(self.weight) <= 1e-7
            self.kp = self.kp[mask]
            self.weight = self.weight[mask]
            self.eig = self.eig[:, mask, :]
            self.occ = self.occ[:, mask, :]
            self.Nk = len(self.kp)

        self.calculate_gap()

    def read_qexml(self, filename=""):

        if filename == "":
            # find a .xml file
            files = os.listdir()
            for f in files:
                if f.endswith('.xml'):
                    filename = f
                    break

        Ha = load_constant("rydberg") * 2.0
        bohr = load_constant("bohr")

        tree = ET.parse(filename)
        band_structure = tree.getroot().find("output").find("band_structure")
        kpoints = band_structure.findall("ks_energies")
        cell = tree.getroot().find("input").find("atomic_structure").find("cell")
        lc = np.zeros((3, 3))
        for ix in range(3):
            lc[ix, :] = [float(x) for x in cell.find(f"a{ix+1}").text.split()]
        lc = lc * bohr
        a1 = np.linalg.norm(lc[0, :])

        self.Nk = len(kpoints)
        lsda_elem = band_structure.find("lsda")
        if lsda_elem is not None and lsda_elem.text.strip().lower() == "true":
            self.Ns = 2
        else:
            self.Ns = 1
        self.Nb = int(kpoints[0].find('eigenvalues').get("size")) // self.Ns

        kp = np.zeros((self.Nk, 3))
        self.weight = np.zeros(self.Nk)
        self.eig = np.zeros((self.Ns, self.Nk, self.Nb))
        self.occ = np.zeros((self.Ns, self.Nk, self.Nb))

        for ik, ks in enumerate(kpoints):
            kp[ik] = [float(kp1) for kp1 in ks.find('k_point').text.split()]
            self.weight[ik] = float(ks.find('k_point').get("weight"))
            eig1 = np.fromstring(ks.find('eigenvalues').text, sep=' ')
            occ1 = np.fromstring(ks.find('occupations').text, sep=' ')
            self.eig[:, ik, :] = eig1.reshape(self.Ns, self.Nb) * Ha
            self.occ[:, ik, :] = occ1.reshape(self.Ns, self.Nb)

        self.kp = kp @ lc.T / a1

        self.calculate_gap()

    def read_wan(self, Nb_pad=0):
        # only support Ns=1 and semiconductor
        files = os.listdir()

        filename = ""
        for f in files:
            if f.endswith('_band.kpt'):
                filename = f
                break
        with open(filename, "r") as f1:
            line = f1.readlines()
        self.Nk = int(line[0].split()[0])
        self.kp = np.zeros((self.Nk, 3))
        self.weight = np.zeros(self.Nk)
        for ik in range(self.Nk):
            word = line[ik+1].split()
            self.kp[ik] = [float(word[0]), float(word[1]), float(word[2])]
            self.weight[ik] = float(word[3])

        filename = ""
        # find a *_band.dat file
        for f in files:
            if f.endswith('_band.dat'):
                filename = f
                break
        with open(filename, "r") as f2:
            line = f2.readlines()

        self.Ns = 1
        energies = []
        current_band = []
        for l in line:
            word = l.split()
            if len(word) == 2:
                current_band.append(float(word[1]))
            elif len(word) == 0 and current_band:
                energies.append(current_band)
                current_band = []
        if current_band:
            energies.append(current_band)

        self.Nb = len(energies)
        self.eig = np.array(energies).T[np.newaxis, :, :]  # (Nb, Nk) -> (Nk, Nb) -> (1, Nk, Nb)
        self.occ = np.zeros((self.Ns, self.Nk, self.Nb))

    def read_parsec(self, lc=None):
        if lc is None:
            lc = np.eye(3)

        bohr = load_constant("bohr")
        pi = load_constant("pi")
        rydberg = load_constant("rydberg")

        if os.path.isfile("bands.dat"):
            filename = "bands.dat"
        elif os.path.isfile("eigen.dat"):
            self.read_parsec_eigen()
            return
        else:
            print("Error: Reading parsec eigenval needs bands.dat or eigen.dat")
            sys.exit()

        with open(filename, "r") as f1:
            line = f1.readlines()

        word = line[1].split()
        self.Ns = int(word[0])
        self.Nb = int(word[1])
        Np = int(word[2])
        ef = float(word[3]) * rydberg
        self.Nk = 0
        for ip in range(Np):
            word = line[ip+3].split()
            self.Nk += int(word[1])

        self.eig = np.zeros((self.Ns, self.Nk, self.Nb))
        self.weight = np.ones(self.Nk)
        kp = np.zeros((self.Nk, 3))

        il = 4 + Np
        for ispin in range(self.Ns):
            for ik in range(self.Nk):
                if ispin == 0:
                    word = line[il].split()
                    kp[ik] = [float(word[3]), float(word[4]), float(word[5])]
                il += 1
                for ib in range(self.Nb):
                    self.eig[ispin, ik, ib] = float(line[il].split()[0])
                    il += 1
        self.eig *= rydberg

        # convert kp from 1/bohr to 1/A
        self.kp = kp @ lc.T / bohr * 0.5 / pi

        self.eig -= ef
        Nvb = np.zeros((self.Ns, self.Nk), dtype=int)
        for ispin in range(self.Ns):
            for ik in range(self.Nk):
                pos_idx = np.where(self.eig[ispin, ik, :] > 1e-6)[0]
                if len(pos_idx) > 0:
                    Nvb[ispin, ik] = pos_idx[0]
                else:
                    Nvb[ispin, ik] = self.Nb

        if (Nvb.max(axis=1) == Nvb.min(axis=1)).all():
            self.is_semic = True
            self.Nvb = Nvb.max(axis=1).tolist()
            # Creates occ if a gap is detected
            self.occ = np.zeros((self.Ns, self.Nk, self.Nb))
            for ispin in range(self.Ns):
                self.occ[ispin, :, 0:self.Nvb[ispin]] = 1.0
        else:
            self.is_semic = False
            self.occ = np.zeros((self.Ns, self.Nk, self.Nb))

        self.calculate_gap()

    def read_parsec_eigen(self):
        # Read the eigenvalues from eigen.dat.
        # kp is missing in this file, so we assume the first k is gamma and read only this point.
        rydberg = load_constant("rydberg")

        filename = "eigen.dat"

        with open(filename, "r") as f0:
            line = f0.readlines()

        self.Nk = 1
        self.Ns = 0
        eig_list = [[], []]
        occ_list = [[], []]
        for l in line:
            word = l.split()
            if len(word) == 0 or word[0][0] in ["!", "#"]:
                continue
            if not word[0].isdigit():
                continue
            if self.Ns == 0:
                if len(word) == 7:
                    self.Ns = 1
                elif len(word) == 8:
                    self.Ns = 2
                else:
                    print("Error: number of columns in eigen.dat should be 7 or 8")
                    sys.exit()
            ik = int(word[5])-1
            if ik > 0:
                continue
            ib = int(word[0])-1
            eig0 = float(word[1])*rydberg
            occ0 = float(word[3])
            if self.Ns == 2 and word[7] == "dn":
                eig_list[1].append(eig0)
                occ_list[1].append(occ0)
            else:  # Ns=1 or spin=up
                eig_list[0].append(eig0)
                occ_list[0].append(occ0)

        self.Nb = len(eig_list[0])
        # (Ns, Nb) -> (Ns, 1, Nb) where Nk=1
        data_eig = eig_list[:self.Ns]
        data_occ = occ_list[:self.Ns]
        self.eig = np.array(data_eig)[:, np.newaxis, :]
        self.occ = np.array(data_occ)[:, np.newaxis, :]

        self.kp = np.zeros((1, 3))
        self.weight = np.ones(1)
        self.calculate_gap()

#########################################################################

    def eig_x(self, kp, rlc=None):
        # transfer the fractional to cartesian in k space
        # if kp is already in cartisian, simply set rlc=I
        if rlc is None:
            rlc = np.eye(3)
        kpc = self.kp @ rlc.T
        kplabelold = ""
        kpout = []
        for ik in range(self.Nk):
            kplabel = kp.findlabel(self.kp[ik], dim=0)
            if kplabel != "elsewhere" and kplabel != kplabelold and (
                kplabelold != "elsewhere" or np.max(np.abs(self.kp[ik] - self.kp[ik-1])) > 0.1
            ):
                kpout.append([0.0])
            else:
                dkpc = float(np.linalg.norm(kpc[ik] - kpc[ik-1]))
                kpout[-1].append(kpout[-1][-1] + dkpc)
            kplabelold = kplabel
        return kpout

    def writegap(self, kp):
        with open("gap.txt", "w") as f0:
            if self.is_semic:
                if self.vbm_k != self.cbm_k:
                    vbm_kl = kp.findlabel(self.kp[self.vbm_k], dim=1)
                    cbm_kl = kp.findlabel(self.kp[self.cbm_k], dim=1)
                    eindg_print = f"Indirect: Eg = {self.eindg:.4f} eV, between {vbm_kl} -> {cbm_kl}"
                    f0.write(eindg_print + "\n")
                    print(eindg_print)
                edg_kl = kp.findlabel(self.kp[self.edg_k], dim=1)
                edg_print = f"Direct:   Eg = {self.edg:.4f} eV, at {edg_kl}"
                f0.write(edg_print + "\n")
                print(edg_print)
            else:
                f0.write("Eg = 0, no band gap\n")

    def calculate_gap(self):
        if self.occ is None or self.eig is None or np.all(self.occ == 0):
            self.is_semic = False
            return

        self.is_semic = not bool(np.any((self.occ > 0.001) & (self.occ < 0.999)))
        if not self.is_semic:
            return

        self.Nvb = []
        for ispin in range(self.Ns):
            unocc = np.where(self.occ[ispin, 0, :] < 0.5)[0]
            if len(unocc) > 0:
                self.Nvb.append(int(unocc[0]))
            else:
                self.Nvb.append(self.Nb)

        self.vbm = -1e10
        self.vbm_k = 0
        self.vbm_s = 0
        self.cbm = 1e10
        self.cbm_k = 0
        self.cbm_s = 0
        # edg: direct band gap
        self.edg = 1e10
        self.edg_k = 0
        self.edg_s = 0

        for ispin in range(self.Ns):
            nvb = self.Nvb[ispin]
            vb_e = self.eig[ispin, :, nvb - 1]
            cb_e = self.eig[ispin, :, nvb]
            dir_gap = cb_e - vb_e

            max_vb_k = int(np.argmax(vb_e))
            if vb_e[max_vb_k] > self.vbm:
                self.vbm = float(vb_e[max_vb_k])
                self.vbm_k = max_vb_k
                self.vbm_s = ispin

            min_cb_k = int(np.argmin(cb_e))
            if cb_e[min_cb_k] < self.cbm:
                self.cbm = float(cb_e[min_cb_k])
                self.cbm_k = min_cb_k
                self.cbm_s = ispin

            min_dir_k = int(np.argmin(dir_gap))
            if dir_gap[min_dir_k] < self.edg:
                self.edg = float(dir_gap[min_dir_k])
                self.edg_k = min_dir_k
                self.edg_s = ispin

        self.eindg = self.cbm - self.vbm
