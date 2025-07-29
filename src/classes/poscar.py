import sys
import os
import numpy as np
import math
import xml.etree.ElementTree as ET
from load_data import load_constant, load_atom_mass


class Poscar:
    """
    This class handles I/O and manipulation of crsytal structures.

    Attributes:
        title (str): system description.
        Ndim (int): dimensionality (3=bulk, 2=slab, 0=molecule).
        lc[3][3] (list of float): lattice vectors (in Angstroms).
        Natom (int): total number of atoms.
        Ntype (int): number of atom types.
        atomtype[Ntype] (list of str): atom type symbols.
        Naint[Ntype] (list of int): number of atoms per type.
        ap[Natom][3] (list of float): atomic positions in fractional coordinates.
        seldyn[Natom][3] (list of bool): optional selective dynamics flags; None if not present.
    """

    def __init__(self):
        self.title = "SYSTEM"
        self.Ndim = 3
        self.lc = []
        self.Natom = 0
        self.Ntype = 0
        self.atomtype = []
        self.Naint = []
        self.ap = []
        self.seldyn = None

    def __str__(self):
        str_out = "POSCAR:\n"
        str_out += " Ndim = " + str(self.Ndim) + "\n"
        str_out += " Natom = " + str(self.Natom) + "\n"
        str_out += " Ntype = " + str(self.Ntype) + "\n"
        for it in range(self.Ntype):
            str_out += "  " + self.atomtype[it] + " " + str(self.Naint[it]) + "\n"
        return str_out

#########################################################################

    def read_vasp(self, filename=""):

        if filename == "":
            filename = "POSCAR"

        f0 = open(filename, "r")
        line = f0.readlines()
        f0.close()

        lineoff = 0
        self.title = line[0].rstrip('\n')
        factor = float(line[1].split()[0])
        for i in range(2, 5):
            word = line[i].split()
            self.lc.append([float(word[0])*factor, float(word[1])*factor, float(word[2])*factor])
        word = line[5].split()
        self.Ntype = len(word)
        for i in range(self.Ntype):
            self.atomtype.append(word[i])
        word = line[6].split()
        for i in range(self.Ntype):
            self.Naint.append(int(word[i]))
            self.Natom = self.Natom+self.Naint[i]
        word = line[7].split()
        if word[0][0].lower() == 's':
            self.seldyn = []
            lineoff = lineoff+1
        word = line[7+lineoff].split()
        if word[0][0].lower() != 'd':
            print("Only atomic position 'Direct' supported")
            sys.exit()
        for i in range(self.Natom):
            word = line[i+lineoff+8].split()
            self.ap.append([float(word[0]), float(word[1]), float(word[2])])
            if self.seldyn is not None:
                self.seldyn.append([bool(word[3]), bool(word[4]), bool(word[5])])
        self.movetobox()
        del line
        del word

    def read_qe(self, filename=""):

        # find a scf.in/nscf.in/relax.in file
        if filename == "":
            files = os.listdir()
            if "scf.in" in files:
                filename = "scf.in"
            elif "nscf.in" in files:
                filename = "nscf.in"
            elif "relax.in" in files:
                filename = "relax.in"
            else:
                print("Error: In read_qe, use input files scf.in, nscf.in, or relax.in, or specify the input file name.")
                sys.exit()

        bohr = load_constant("bohr")

        f1 = open(filename, "r")
        line = f1.readlines()
        f1.close()

        Fat = False

        for l in range(len(line)):
            word = line[l].split()
            if len(word) >= 2 and word[0] == "CELL_PARAMETERS":
                word[1] = word[1].replace("(", "").replace("{", "")
                if word[1][0].lower() == "a":
                    factor = 1.0
                elif word[1][0].lower() == "b":
                    factor = bohr
                else:
                    print(word[1])
                    print("Only Angstrom and Bohr are supported.")
                    sys.exit()
                for il in range(3):
                    wordlc = line[l+il+1].split()
                    self.lc.append([float(wordlc[0]), float(wordlc[1]), float(wordlc[2])])
                continue
            elif len(word) >= 2 and word[0] == "ATOMIC_POSITIONS":
                Fat = True
                continue

            if Fat == True:
                try:
                    self.ap.append([float(word[1]), float(word[2]), float(word[3])])
                except:
                    Fat = False
                    continue

                self.Natom += 1
                if self.Ntype == 0 or self.atomtype[-1] != word[0]:
                    self.Ntype += 1
                    self.atomtype.append(word[0])
                    self.Naint.append(1)
                else:
                    self.Naint[-1] += 1

        for i in range(3):
            for j in range(3):
                self.lc[i][j] = self.lc[i][j]*factor

        self.movetobox()

    def read_xml(self, filename=""):

        if filename == "":
            # find a .xml file
            files = os.listdir()
            for f in files:
                if f.endswith('.xml'):
                    filename = f
                    break
        if filename == "":
            print("Error: .xml file is not found")
            sys.exit()

        bohr = load_constant("bohr")

        tree = ET.parse(filename)
        cell = tree.getroot().find("output").find("atomic_structure").find("cell")

        for ix in range(3):
            self.lc.append([float(lc1) for lc1 in cell.find("a"+str(ix+1)).text.split()])
        atoms = tree.getroot().find("output").find("atomic_structure").find("atomic_positions").findall("atom")

        for atom in atoms:
            atom_name = atom.get("name")
            if len(self.atomtype) == 0 or self.atomtype[-1] != atom_name:
                self.atomtype.append(atom_name)
                self.Naint.append(1)
            else:
                self.Naint[-1] += 1
            self.ap.append([float(ap) for ap in atom.text.split()])

        # convert Cartesian to direct
        self.ap = (np.array(self.ap)@np.linalg.inv(np.array(self.lc))).tolist()
        for ix in range(3):
            for iy in range(3):
                self.lc[ix][iy] = self.lc[ix][iy]*bohr

        self.Natom = len(self.ap)
        self.Ntype = len(self.atomtype)

    def read_prt(self, filename=""):

        if filename == "":
            filename = "input"

        bohr = load_constant("bohr")

        f1 = open(filename, "r")
        line = f1.readlines()
        f1.close()

        Flc = False
        Fat = False

        for l in range(len(line)):
            word = line[l].split()
            if len(word) >= 2 and word[0] == "begin" and word[1] == "latticevecs":
                Flc = True
                continue
            elif len(word) >= 2 and word[0] == "end" and word[1] == "latticevecs":
                Flc = False
                continue
            elif len(word) >= 2 and word[0] == "begin" and word[1] == "coordinates":
                Fat = True
                continue
            elif len(word) >= 2 and word[0] == "end" and word[1] == "coordinates":
                Fat = False
                continue

            if Flc == True:
                if len(word) >= 4 and word[0] == "coord":
                    self.lc.append([float(word[1]), float(word[2]), float(word[3])])
                elif len(word) >= 2 and word[0] == "volume":
                    volume = float(word[1])
            elif Fat == True:
                if len(word) >= 2 and word[0] == "newtype":
                    self.Ntype += 1
                    self.atomtype.append(word[1])
                    self.Naint.append(0)
                elif len(word) >= 4 and word[0] == "coord":
                    self.ap.append([float(word[1]), float(word[2]), float(word[3])])
                    self.Naint[-1] += 1
                    self.Natom += 1

        detlc = self.volume()
        factor = (volume/detlc)**(1.0/3.0)*bohr

        for i in range(3):
            for j in range(3):
                self.lc[i][j] = self.lc[i][j]*factor

        self.movetobox()

    def read_parsec(self, filename=""):

        if filename == "":
            filename = "parsec.in"

        bohr = load_constant("bohr")

        f1 = open(filename, "r")
        line = f1.readlines()
        f1.close()

        self.Ndim = 0
        for il in range(len(line)):
            word = line[il].replace(":", " ").replace("=", " ").split()
            if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                continue
            if word[0].lower() == "boundary_conditions":
                if word[1].lower() == "slab":
                    self.Ndim = 2
                    break
                if word[1].lower() == "bulk":
                    self.Ndim = 3
                    break

        Fat = False
        Flc = False
        self.lc = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        factor_lc = bohr
        for il in range(len(line)):
            word = line[il].replace(":", " ").replace("=", " ").split()
            if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                continue
            if len(word) >= 2 and word[0].lower() == "begin" and word[1].lower() == "atom_coord":
                Fat = True
                rlc = np.linalg.inv(np.array(self.lc)).tolist()
                continue
            if len(word) >= 2 and word[0].lower() == "end" and word[1].lower() == "atom_coord":
                Fat = False
                continue
            if len(word) >= 2 and word[0].lower() == "begin" and word[1].lower() == "cell_shape":
                Flc = True
                ilc = 0
                continue
            if len(word) >= 2 and word[0].lower() == "end" and word[1].lower() == "cell_shape":
                Flc = False
                continue
            if word[0].lower() == "boundary_sphere_radius":
                alat = float(word[1]) * 2.0  # in Bohr here
                if len(word) < 3 or word[2].lower() != "ang":
                    alat = alat * bohr
                if self.Ndim == 0:
                    self.lc[0][0] = alat
                    self.lc[1][1] = alat
                    self.lc[2][2] = alat
                elif self.Ndim == 2:
                    self.lc[2][2] = alat
            if word[0].lower() == "lattice_vector_scale":
                factor_lc = float(word[1])
                if len(word) < 3 or word[2].lower() != "ang":
                    factor_lc = factor_lc * bohr
            if word[0].lower() == "atom_types_num":
                self.Ntype = int(word[1])
            if word[0].lower() == "atom_type":
                self.atomtype.append(word[1])
                self.Naint.append(0)
            if word[0].lower() == "coordinate_unit":
                if word[1].lower() == "cartesian_ang":
                    lcart = True
                    factor_ap = 1.0
                elif word[1].lower() == "lattice_vectors":
                    lcart = False
                    factor_ap = 1.0
                else:  # cartesian_bohr
                    lcart = True
                    factor_ap = bohr
            if Flc:
                for ix in range(3):
                    self.lc[ilc][ix] = float(word[ix]) * factor_lc
                ilc += 1
            if Fat:
                ap = []
                for ix in range(3):
                    ap.append(float(word[ix])*factor_ap)
                self.ap.append(ap)
                self.Naint[-1] += 1
                self.Natom += 1
        if lcart:
            self.ap = (np.array(self.ap)@np.linalg.inv(np.array(self.lc))).tolist()
        if self.Ndim == 0:
            shift = [0.5, 0.5, 0.5]
        elif self.Ndim == 2:
            shift = [0.0, 0.0, 0.5]
        else:  # bulk
            shift = [0.0, 0.0, 0.0]
        for ia in range(self.Natom):
            for ix in range(3):
                self.ap[ia][ix] += shift[ix]

    def read_xyz(self, filename=""):

        if filename == "":
            # find a .xyz file
            files = os.listdir()
            for f in files:
                if f.endswith('.xyz'):
                    filename = f
                    break
        if filename == "":
            print("Error: .xyz file is not found")
            sys.exit()

        bohr = load_constant("bohr")

        self.Ndim = 0
        f0 = open(filename, "r")
        line = f0.readlines()
        f0.close()

        self.Natom = int(line[0].split()[0])
        lbohr = False
        word = line[1].split()
        ap_list = []
        if len(word) >= 1 and word[0][0] in ["B", "b"]:
            lbohr = True
        if len(line) < self.Natom+2:
            print("Error: Length of file "+filename+" is less than Natom+2 "+str(self.Natom))
            sys.exit()
        for il in range(2, self.Natom+2):
            word = line[il].split()
            if self.Ntype == 0 or self.atomtype[-1] != word[0]:
                self.Ntype += 1
                self.atomtype.append(word[0])
                self.Naint.append(1)
            else:
                self.Naint[-1] += 1
            ap_list.append([float(word[1]), float(word[2]), float(word[3])])

        ap = np.array(ap_list)
        del ap_list

        if lbohr:
            ap = ap * bohr

        radius = np.max(np.linalg.norm(ap, axis=1))
        lc0 = 2.0 * (radius+5.0)  # vacuum is set to 5 A
        self.lc = [[lc0, 0.0, 0.0], [0.0, lc0, 0.0], [0.0, 0.0, lc0]]

        self.ap = (ap / lc0 + 0.5).tolist()
        del ap

########################################################################

    def write_vasp(self, filename="POSCAR.new"):
        f1 = open(filename, "w")
        f1.write(self.title+'\n')
        f1.write('1.0\n')
        for i in range(3):
            for j in range(3):
                f1.write(f"{self.lc[i][j]:21.12f}")
            f1.write('\n')
        for i in range(self.Ntype):
            f1.write(self.atomtype[i]+"   ")
        f1.write('\n')
        for i in range(self.Ntype):
            f1.write(str(self.Naint[i])+"  ")
        f1.write('\n')
        if self.seldyn is not None:
            f1.write('Selective dynamics\n')
        f1.write('Direct\n')
        for i in range(self.Natom):
            for j in range(3):
                f1.write(f"{self.ap[i][j]:18.12f}")
            if self.seldyn is not None:
                for j in range(3):
                    if self.seldyn[i][j]:
                        f1.write(" T")
                    else:
                        f1.write(" F")
            f1.write('\n')
        f1.close()

    def write_qe(self, filename="qe_st.dat"):
        mass = load_atom_mass()

        f2 = open(filename, "w")
        f2.write("CELL_PARAMETERS angstrom\n")
        for ix1 in range(3):
            for ix2 in range(3):
                f2.write(f"{self.lc[ix1][ix2]:18.12f}")
            f2.write("\n")

        f2.write("ATOMIC_SPECIES\n")
        for i in range(self.Ntype):
            f2.write("  " + self.atomtype[i] + "  " + str(mass[self.atomtype[i]]) + "  " + self.atomtype[i] + ".upf\n")
        f2.write("ATOMIC_POSITIONS crystal\n")
        ij = 0
        ik = 0
        for i in range(self.Natom):
            f2.write(f"  {self.atomtype[ij]:2s}{self.ap[i][0]:18.12f}{self.ap[i][1]:18.12f}{self.ap[i][2]:18.12f}\n")
            ik = ik + 1
            if ik == self.Naint[ij]:
                ij = ij + 1
                ik = 0

        k_grid = self.k_grid()
        f2.write("K_POINTS automatic\n")
        f2.write(f"  {k_grid[0]:d}  {k_grid[1]:d}  {k_grid[2]:d} 0 0 0\n")
        f2.close()

    def write_prt(self, filename="prt_st.dat"):
        bohr = load_constant("bohr")

        volume = self.volume()*bohr**(-3)
        # convert Angstrom^3 to bohr^3

        f2 = open(filename, "w")
        f2.write("begin latticevecs\n")
        for i in range(3):
            f2.write(f"coord{self.lc[i][0]:18.12f}{self.lc[i][1]:18.12f}{self.lc[i][2]:18.12f}\n")
        f2.write(f"volume {volume:24.16f}\nend latticevecs\n\nbegin coordinates\n")
        k = 0
        for i in range(self.Ntype):
            f2.write("newtype "+self.atomtype[i]+"\n")
            for j in range(self.Naint[i]):
                f2.write(f"coord{self.ap[k][0]:18.12f}{self.ap[k][1]:18.12f}{self.ap[k][2]:18.12f}\n")
                k = k+1
        f2.write("end coordinates\n\n")
        f2.close()

    def write_parsec(self, filename="parsec_st.dat", lcartesian=False, lbohr=False):
        bohr = load_constant("bohr")

        if self.Ndim < 3:
            lcartesian = True
        f2 = open(filename, "w")
        if lbohr:
            funit = bohr
        else:
            funit = 1.0

        if self.Ndim == 0:
            radius = max(self.lc[0][0], self.lc[1][1], self.lc[2][2])*0.5/funit
            if lbohr:
                f2.write(f"boundary_sphere_radius {radius:.12g}\n\n")
            else:
                f2.write(f"boundary_sphere_radius {radius:.12g} ang\n\n")
        elif self.Ndim == 1:
            f2.write("boundary_conditions wire\n")
            radius = max(self.lc[1][1], self.lc[2][2])*0.5/funit
            if lbohr:
                f2.write(f"boundary_sphere_radius {radius:.12g}\n\n")
            else:
                f2.write(f"boundary_sphere_radius {radius:.12g} ang\n\n")
        elif self.Ndim == 2:
            f2.write("boundary_conditions slab\n")
            radius = self.lc[2][2]*0.5/funit
            if lbohr:
                f2.write(f"boundary_sphere_radius {radius:.12g}\n\n")
            else:
                f2.write(f"boundary_sphere_radius {radius:.12g} ang\n\n")
        else:  # bulk
            f2.write("boundary_conditions bulk\n")

        if self.Ndim > 0:
            if not lbohr:
                f2.write("lattice_vector_scale 1.0 ang\n")
            f2.write("begin cell_shape\n")
            for ix1 in range(self.Ndim):
                for ix2 in range(3):
                    f2.write(f"{self.lc[ix1][ix2]/funit:18.12f}")
                f2.write("\n")
            f2.write("end cell_shape\n\n")

            f2.write("kpoint_method mp\n\n")
            f2.write("begin monkhorst_pack_grid\n")
            k_grid = self.k_grid()
            for ix in range(self.Ndim):
                f2.write(f"  {k_grid[ix]:d}")
            f2.write("\nend monkhorst_pack_grid\n\n")
            f2.write("begin monkhorst_pack_shift\n")
            f2.write("0.0  0.0  0.0\n")
            f2.write("end monkhorst_pack_shift\n\n")

        f2.write("atom_types_num " + str(len(self.atomtype)) + "\n")
        if lcartesian:
            if lbohr:
                f2.write("coordinate_unit cartesian_bohr\n\n")
            else:
                f2.write("coordinate_unit cartesian_ang\n\n")
            apc = self.cartesian(factor=1.0/funit)
        else:
            f2.write("coordinate_unit lattice_vectors\n\n")
        k = 0
        for i in range(self.Ntype):
            f2.write("atom_type " + self.atomtype[i] + "\n")
            # f2.write("Pseudopotential_Format: \n")
            # f2.write("Core_Cutoff_Radius: \n")
            f2.write("local_component s\n")  # to be modified
            # f2.write("Potential_Num: \n\n")
            # f2.write("begin Electron_Per_Orbital\n")
            # f2.write("# S P D F\n \n")
            # f2.write("end Electron_Per_Orbital\n\n")

            f2.write("begin atom_coord\n")
            if lcartesian:
                for ia in range(self.Naint[i]):
                    f2.write(f"{apc[k][0]:18.12f}{apc[k][1]:18.12f}{apc[k][2]:18.12f}\n")
                    k = k + 1
            else:
                for ia in range(self.Naint[i]):
                    f2.write(f"{self.ap[k][0]:18.12f}{self.ap[k][1]:18.12f}{self.ap[k][2]:18.12f}\n")
                    k = k + 1
            f2.write("end atom_coord\n\n")

        f2.close()

    def write_wannier90(self, filename="wannier90_st.dat"):
        f2 = open(filename, "w")
        f2.write("Begin Unit_Cell_Cart\n")
        for i in range(3):
            f2.write(f"{self.lc[i][0]:18.12f}{self.lc[i][1]:18.12f}{self.lc[i][2]:18.12f}\n")
        f2.write("End Unit_Cell_Cart\n\n")

        f2.write("Begin Projections\n")
        f2.write("  random\n")
        for i in range(self.Ntype):
            f2.write("  "+self.atomtype[i]+"  :\n")
        f2.write("End Projections\n\n")

        f2.write("Begin Atoms_Frac\n")
        ij = 0
        ik = 0
        for i in range(self.Natom):
            f2.write(f"  {self.atomtype[ij]:2s}{self.ap[i][0]:18.12f}{self.ap[i][1]:18.12f}{self.ap[i][2]:18.12f}\n")
            ik = ik+1
            if ik == self.Naint[ij]:
                ij = ij+1
                ik = 0
        f2.write("End Atoms_Frac\n\n")
        f2.close()

    def write_xyz(self, filename="structure.xyz"):

        self.Ndim = 0

        f2 = open(filename, "w")
        f2.write(str(self.Natom) + "\n\n")
        it1 = 0
        it2 = 0
        apc = self.cartesian()
        for ia in range(self.Natom):
            # write only 10 digits after the decimal point to eliminate small residuals
            f2.write(f"  {self.atomtype[it1]:2s}{apc[ia][0]:18.12f}{apc[ia][1]:18.12f}{apc[ia][2]:18.12f}\n")
            it2 = it2 + 1
            if it2 == self.Naint[it1]:
                it1 = it1 + 1
                it2 = 0

        f2.close()

########################################################################

    def movetobox(self):
        tol = 1e-8
        ap = np.array(self.ap)
        ap = ap-ap//1.0
        ap[(ap < tol) | (ap > 1.0-tol)] = 0.0
        self.ap = ap.tolist()

    def move(self, disp, lcart=False, lbox=True):
        # move atoms by a set of vectors
        # disp[Na][3] or disp[3]
        # if lcart=True, the vector is defined in a cartesian coordinate system using angstrom
        # if lcart=False, the vector is defined in the lattice coordinate
        # if lbox=True, the atoms will be moved into the box

        disp = np.array(disp)
        shape = disp.shape
        if not (shape == (3,) or shape == (self.Natom, 3)):
            print("Error: the dimension of displacement vector should be (3) or (Na,3)")
            sys.exit()
        if lcart:
            # if the disp is in cartesian, convert it to direct
            rlc = self.rlc()
            disp = np.array(disp)@np.array(rlc)
        self.ap = (np.array(self.ap)+disp).tolist()
        if lbox:
            self.movetobox()

    def rotate(self, theta):
        # rotate the structure around the z-axis (in degree)
        t = math.radians(theta)
        c = math.cos(t)
        s = math.sin(t)
        M = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
        self.lc = (np.array(self.lc)@M).tolist()

    def flip(self):
        # flip the z coordinate of the structure
        ap = np.array(self.ap)
        ap[:, 2] = 1.0-ap[:, 2]
        self.ap = ap.tolist()
        self.movetobox()

    def supercell(self, N):
        # N should be a list of 3 integers
        N_total = N[0]*N[1]*N[2]

        M = np.array([[N[0], 0, 0], [0, N[1], 0], [0, 0, N[2]]])
        self.lc = (M@np.array(self.lc)).tolist()

        shift = np.zeros((N_total, 3))
        i = 0
        for iz in range(N[2]):
            for iy in range(N[1]):
                for ix in range(N[0]):
                    shift[i, :] = np.array([ix, iy, iz])/np.array(N)
                    i += 1

        ap = np.array(self.ap)
        apnew = np.zeros((N_total*self.Natom, 3))
        for ia in range(self.Natom):
            apnew[ia*N_total:(ia+1)*N_total, :] = ap[ia, :]/np.array(N)+shift
        self.ap = apnew.tolist()

        for itype in range(self.Ntype):
            self.Naint[itype] = self.Naint[itype]*N_total
        self.Natom = self.Natom*N_total

        self.movetobox()

    def vacuum(self, z_vac):
        # Add a vacuum layer to the structure
        # z_vac is in angstrom
        # This version assumes a1 x a2 // a3 // z
        if z_vac <= 0:
            print("Error: z_vac <= 0")
            sys.exit()

        a3_new = self.lc[2][2] + z_vac

        self.movetobox()

        ia = 0
        ap_new = []
        for itype in range(self.Ntype):
            ap_new0 = []
            for iaint in range(self.Naint[itype]):
                if self.ap[ia][2] < 1e-8:
                    ap_new0.append([self.ap[ia][0], self.ap[ia][1], 1.0])
                ia += 1
            ap_new.append(ap_new0)

        for itype in range(self.Ntype):
            for ap in ap_new[itype]:
                self.add_atom(itype, ap)

        for ia in range(self.Natom):
            self.ap[ia][2] = (self.ap[ia][2] - 0.5) * self.lc[2][2] / a3_new + 0.5

        self.lc[2][2] = a3_new

    def add_atom(self, itype, ap, new_type=None, add_to_head=False):
        # If new_type is defined, a new atom type will be inserted.
        # add_to_head defines if the atom will be added to the head.
        if itype < 0:
            print("Error: itype ("+str(itype)+") < 0 in add_atom function.")
            sys.exit()
        elif itype > self.Ntype:
            print("Error: itype ("+str(itype)+") > self.Ntype ("+str(self.Ntype)+") in add_atom function.")
            sys.exit()

        if new_type is None:
            if itype == self.Ntype:
                print("Warning: itype ("+str(itype)+") > self.Ntype ("+str(self.Ntype)+") in add_atom function,")
                print("         but atom name is not defined.")
                new_type = "X"
            else:  # This is not a new type
                self.Natom += 1
                if add_to_head:
                    self.ap.insert(sum(self.Naint[:itype]), ap)
                else:
                    self.ap.insert(sum(self.Naint[:itype+1]), ap)
                self.Naint[itype] += 1
                return

        else:  # This is a new type
            self.Natom += 1
            self.Ntype += 1
            self.atomtype.insert(itype, new_type)
            self.Naint.insert(itype, 1)
            self.ap.insert(sum(self.Naint[:itype]), ap)

    def delete_atom(self, ia):
        if ia < 0 or ia >= self.Natom:
            print("Error: ia ("+ia+") is out of range")
            sys.exit()

        self.Natom -= 1
        del self.ap[ia]
        iat = 0
        for itype in range(self.Ntype):
            iat = iat + self.Naint[itype]
            if iat > ia:
                self.Naint[itype] -= 1
                if self.Naint[itype] == 0:
                    del self.Naint[itype]
                    del self.atomtype[itype]
                    self.Ntype -= 1
                break

    def replace_atom(self, ia, new_type, add_to_head=True):
        if ia < 0 or ia >= self.Natom:
            print("Error: ia ("+ia+") is out of range")
            sys.exit()

        ap = self.ap[ia]
        self.delete_atom(ia)

        for itype in range(self.Ntype):
            if self.atomtype[itype] == new_type:
                self.add_atom(itype, ap, add_to_head=add_to_head)
                return
        if add_to_head:
            self.add_atom(0, ap, new_type)
        else:
            self.add_atom(self.Ntype, ap, new_type)

    def new_lc(self, lc_new):
        apc = np.array(self.cartesian())
        self.ap = (apc @ np.linalg.inv(lc_new)).tolist()
        self.lc = lc_new

    def rlc(self):
        pi = load_constant("pi")
        return (np.linalg.inv(self.lc)*(2.0*pi)).tolist()

    def volume(self):
        return abs(np.linalg.det(self.lc))

    def cartesian(self, shift=None, factor=1.0):
        if shift is None:
            if self.Ndim == 3:
                shift = [0.0, 0.0, 0.0]
            elif self.Ndim == 2:
                shift = [0.0, 0.0, -0.5]
            elif self.Ndim == 0:
                shift = [-0.5, -0.5, -0.5]
        return ((np.array(self.ap) + np.array(shift)) @ np.array(self.lc) * factor).tolist()

    def atom_list(self):
        atom = []
        for itype in range(self.Ntype):
            for iaint in range(self.Naint[itype]):
                atom.append(self.atomtype[itype])
        return atom

    def read_atom_list(self, atom_string):
        atom_flag = np.zeros(self.Natom, dtype=bool)
        if atom_string == "all":
            atom_flag[:] = True
        else:
            for atom_section in atom_string.split(","):
                is_atom_name = False
                ia_start = 0
                ia_end = 0
                for itype in range(len(self.atomtype)):
                    ia_end = ia_end + self.Naint[itype]
                    if atom_section == self.atomtype[itype]:
                        atom_flag[ia_start: ia_end] = True
                        is_atom_name = True
                    ia_start = ia_end
                if is_atom_name == False:
                    atom_list = atom_section.replace("~", "-").split("-")
                    if len(atom_list) == 1:
                        atom_flag[int(atom_list[0])-1] = True
                    else:
                        atom_flag[int(atom_list[0]) - 1: int(atom_list[1])] = True
        return atom_flag

    def k_grid(self):
        grid = [1]*self.Ndim
        for ix in range(self.Ndim):
            grid[ix] = math.ceil(30.0/(np.linalg.norm(self.lc[ix])))
        return grid
