import sys
import os
import re
import numpy as np
import xml.etree.ElementTree as ET
from load_data import load_constant


class Eigenval:
    """
    Parses and analyzes band structure eigenvalue data and orbital projections.

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
        Na (int): The number of atoms.
        Norb (int): The number of orbitals.
        proj[Ns, Nk, Nb, Na, Norb] (numpy.ndarray): Orbital projections.
        complex[Ns, Nk, Nb, Na, Norb] (numpy.ndarray): Complex projections with phase.
        orb_name[Norb] (list of str): Names of the orbitals.
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
        # Projection attributes
        self.Na = 0
        self.Norb = 0
        self.proj = None
        self.complex = None
        self.orb_name = []

    def __str__(self):
        str_out = "EIGENVAL:\n"
        str_out += f" Nk = {self.Nk}\n"
        str_out += f" Nb = {self.Nb}\n"
        str_out += f" Ns = {self.Ns}\n"
        if self.proj is not None or self.complex is not None:
            str_out += " Contains projections\n"
            str_out += f" Na = {self.Na}\n"
            str_out += f" Norb = {self.Norb}\n"
        else:
            str_out += " Does not contain projections\n"
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

    def read_vasp_procar(self, filename="PROCAR", update_eig=None):
        with open(filename, "r") as f0:
            line = f0.readlines()

        word = line[0].split()
        has_complex = False
        if len(word) >= 5 and word[4].startswith("phase"):
            has_complex = True

        word = line[1].split()
        Nk = int(word[3])
        Nb = int(word[7])
        self.Na = int(word[11])

        Ns = 1
        if not has_complex and len(line) > (((self.Na + 5) * Nb + 3) * Nk + 1) * 1.5:
            # the file length should be (((Na+5)Nb+3)Nk+1)Ns+1
            Ns = 2
        if has_complex and len(line) > (((2 * self.Na + 7) * Nb + 3) * Nk + 1) * 1.5:
            # the file length should be (((2*Na+7)Nb+3)Nk+1)Ns+1
            Ns = 2

        word = line[7].split()
        self.Norb = len(word) - 2

        print(f"Ns: {Ns}  Na: {self.Na}  Nk: {Nk}  Nb: {Nb}  Norb: {self.Norb}")

        if update_eig is not None:
            should_update_eig = update_eig
        else:
            should_update_eig = self.eig is None

        if should_update_eig:
            self.Ns = Ns
            self.Nk = Nk
            self.Nb = Nb
            self.kp = np.zeros((self.Nk, 3))
            self.weight = np.zeros(self.Nk)
            self.eig = np.zeros((self.Ns, self.Nk, self.Nb))
            self.occ = np.zeros((self.Ns, self.Nk, self.Nb))

        self.proj = np.zeros((Ns, Nk, Nb, self.Na, self.Norb))
        if has_complex:
            self.complex = np.zeros((Ns, Nk, Nb, self.Na, self.Norb), dtype=np.complex128)

        ispin = -1
        for this_line in line:
            word = this_line.split()
            if len(word) == 0 or word[0][0] == "!":
                continue
            if len(word) >= 9 and word[0] == "#" and word[2] == "k-points:":
                ispin += 1
            elif word[0] == "k-point":
                ik = int(word[1]) - 1
                if should_update_eig:
                    self.weight[ik] = float(word[-1])
                    if len(word) >= 6 and word[2] == ":":
                        self.kp[ik] = [float(word[3]), float(word[4]), float(word[5])]
            elif word[0] == "band":
                ib = int(word[1]) - 1
                if should_update_eig:
                    self.eig[ispin, ik, ib] = float(word[4])
                    self.occ[ispin, ik, ib] = float(word[7])
            elif word[0].isdigit():
                ia = int(word[0]) - 1
                if self.complex is not None and len(word) > 2 * self.Norb:
                    for iorb in range(self.Norb):
                        self.complex[ispin, ik, ib, ia, iorb] = complex(
                            float(word[iorb * 2 + 1]),
                            float(word[iorb * 2 + 2])
                        )
                else:
                    for iorb in range(self.Norb):
                        self.proj[ispin, ik, ib, ia, iorb] = float(word[iorb + 1])

        if self.Norb == 9:
            self.orb_name = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "x2-y2"]
        elif self.Norb == 16:
            self.orb_name = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz",
                             "x2-y2", "fy3x2", "fxyz", "fyz2", "fz3", "fxz2", "fzx2", "fx3"]
        else:
            print(f"Norb = {self.Norb} is not support yet")
            sys.exit()

        if should_update_eig:
            self.calculate_gap()

    def read_qe(self, filename=""):

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

    def read_qe_projwfc(self, filename_out=None, has_complex=False, filename_xml=None):

        if filename_out is None:
            filename_out = "projwfc.out"

        if has_complex and filename_xml is None:
            filename_xml = "pwscf.save/atomic_proj.xml"

        with open(filename_out, "r") as f0:
            line = f0.readlines()

        l_map = []
        orb_map = []
        atom_map = []
        is_proj = False
        is_first_k = True
        Ns = 1
        for this_line in line:
            word = this_line.split()
            if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                continue
            if word[0] == "state":
                match = re.search(r"atom\s+(\d+).*l=(\d+)\s+m=\s*(\d+)", this_line)
                ia = int(match.group(1)) - 1
                il = int(match.group(2))
                im = int(match.group(3))
                l_map.append(il)
                orb_map.append(il**2 + im - 1)
                atom_map.append(ia)
            elif word[0] == "natomwfc":
                Nproj = int(word[2])
            elif word[0] == "nkstot":
                nkstot = int(word[2])  # nkstot = Ns * Nk (total number of k-states)
            elif word[0] == "nbnd":
                Nb = int(word[2])
            elif word[0] == "k":
                k_point = np.array(word[2:5], dtype=float)
                if is_first_k:
                    is_first_k = False
                    self.Na = max(atom_map) + 1
                    first_k_point = k_point.copy()
                    self.Norb = (max(l_map) + 1)**2
                    if not has_complex:
                        proj0 = np.zeros((nkstot, Nb, self.Na, self.Norb))
                    ik = -1
                ik += 1
                ib = -1
            elif word[0] == "psi":
                ib += 1
                is_proj = True
            elif word[0] == "|psi|^2":
                is_proj = False
            elif len(word) >= 2 and word[0] == "spin" and word[1] == "down":
                Ns = 2

            if not has_complex and is_proj:
                number = re.findall(r"[-+]?(?:\d*\.*\d+)", this_line)  # find all numbers
                for ii in range(0, len(number), 2):
                    # proj0[k][b][a][orb]
                    proj0[ik, ib, atom_map[int(number[ii + 1]) - 1],
                          orb_map[int(number[ii + 1]) - 1]] += float(number[ii])

        Nk = nkstot // Ns
        print(f"Ns: {Ns}  Na: {self.Na}  Nk: {Nk}  Nb: {Nb}  Norb: {self.Norb}")
        if self.Nk == 0:
            self.Ns = Ns
            self.Nk = Nk
            self.Nb = Nb

        if not has_complex:
            self.proj = proj0.reshape(Ns, Nk, Nb, self.Na, self.Norb)
        else:
            root = ET.parse(filename_xml).getroot()
            all_projs = root.find("EIGENSTATES").findall("PROJS")

            self.complex = np.zeros((Ns, Nk, Nb, self.Na, self.Norb), dtype=np.complex128)
            for igroup, projs in enumerate(all_projs):
                ispin = igroup // Nk
                ik = igroup % Nk
                for wfc in projs:
                    iwfc = int(wfc.get("index")) - 1
                    ia = atom_map[iwfc]
                    iorb = orb_map[iwfc]
                    cvals = np.fromstring(wfc.text, sep=" ").view(np.complex128)
                    self.complex[ispin, ik, :, ia, iorb] = cvals

            self.proj = np.abs(self.complex)**2

        if self.Norb == 4:
            self.orb_name = ["s", "pz", "px", "py"]
        elif self.Norb == 9:
            self.orb_name = ["s", "pz", "px", "py", "dz2", "dxz", "dyz", "x2-y2", "dxy"]
        elif self.Norb == 16:
            self.orb_name = ["s", "pz", "px", "py", "dz2", "dxz", "dyz", "x2-y2",
                             "dxy", "fz3", "fxz2", "fyz2", "fzx2", "fxyz", "fx3", "fy3x2"]
        else:
            print(f"Norb = {self.Norb} is not support yet")
            sys.exit()

    def read_wannier90(self, Nb_pad=0):
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
            word = line[ik + 1].split()
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
            word = line[ip + 3].split()
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
            ik = int(word[5]) - 1
            if ik > 0:
                continue
            ib = int(word[0]) - 1
            eig0 = float(word[1]) * rydberg
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
        kp_label_old = ""
        kp_out = []
        for ik in range(self.Nk):
            kp_label = kp.findlabel(self.kp[ik], dim=0)
            is_new_label = (kp_label != "elsewhere") and (kp_label != kp_label_old)
            is_jump = (kp_label_old != "elsewhere") or (np.max(np.abs(self.kp[ik] - self.kp[ik - 1])) > 0.1)

            if is_new_label and is_jump:
                kp_out.append([0.0])
            else:
                dkpc = float(np.linalg.norm(kpc[ik] - kpc[ik - 1]))
                kp_out[-1].append(kp_out[-1][-1] + dkpc)
            kp_label_old = kp_label
        return kp_out

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

    def plot_proj(self, atom_flag, orb_flag):
        # plot_proj[Ns][Nb][Nk]  # numpy
        if self.proj is None:
            print("Error: proj data not loaded.")
            return None
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
                print(f"projector {orb} does not exist")

        for orb in orb_list:
            if orb in self.orb_name:
                orb_flag[self.orb_name.index(orb)] = True
            else:
                print(f"projector {orb} does not exist")

        return orb_flag

    def calculate_pdos(self, e_pdos, sigma):
        # pdos[Ns, Ne, Na, Norb]
        if self.proj is None:
            print("Error: proj data not loaded.")
            return None

        gaussian_coeff = (1 / (sigma * np.sqrt(2 * np.pi)))
        delta = (e_pdos[None, :, None, None] - self.eig[:, None, :, :]) / sigma  # [Ns, Ne, Nk, Nb]
        smearing = gaussian_coeff * np.exp(-0.5 * delta**2)

        weighted_proj = self.proj * self.weight[None, :, None, None, None]

        pdos = np.einsum('sekb,skbao->seao', smearing, weighted_proj)

        return pdos
