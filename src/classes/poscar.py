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
        lc[3, 3] (numpy.ndarray): lattice vectors (in Angstroms).
        Natom (int): total number of atoms.
        Ntype (int): number of atom types.
        atomtype[Ntype] (list of str): atom type symbols.
        Naint[Ntype] (list of int): number of atoms per type.
        ap[Natom, 3] (numpy.ndarray): atomic positions in fractional coordinates.
        seldyn[Natom][3] (list of bool): optional selective dynamics flags; None if not present.
    """

    def __init__(self):
        self.title = "SYSTEM"
        self.Ndim = 3
        self.lc = np.zeros((3, 3))
        self.Natom = 0
        self.Ntype = 0
        self.atomtype = []
        self.Naint = []
        self.ap = np.empty((0, 3))
        self.seldyn = None
        self.k_grid = None
        self.k_grid_shift = None

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
        for il in range(3):
            word = line[il+2].split()
            for ix in range(3):
                self.lc[il, ix] = float(word[ix])
        self.lc = self.lc * factor

        word = line[5].split()
        self.Ntype = len(word)
        for i in range(self.Ntype):
            self.atomtype.append(word[i])
        word = line[6].split()
        for i in range(self.Ntype):
            self.Naint.append(int(word[i]))
        self.Natom = sum(self.Naint)
        word = line[7].split()
        if word[0][0].lower() == 's':
            self.seldyn = []
            lineoff = lineoff+1
        word = line[7+lineoff].split()
        if word[0][0].lower() != 'd':
            print("Only atomic position 'Direct' supported")
            sys.exit()
        self.ap = np.zeros((self.Natom, 3))
        for ia in range(self.Natom):
            word = line[ia+lineoff+8].split()
            for ix in range(3):
                self.ap[ia, ix] = float(word[ix])
            if self.seldyn is not None:
                self.seldyn.append([bool(word[3]), bool(word[4]), bool(word[5])])
        self.wrap_to_cell()
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
                print("Error: In read_qe, use input files scf.in, nscf.in, or relax.in,\n")
                print("       or specify the input file name.")
                sys.exit()

        bohr = load_constant("bohr")

        f1 = open(filename, "r")
        line = f1.readlines()
        f1.close()

        Fat = False
        ap = []
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
                    for ix in range(3):
                        self.lc[il, ix] = float(wordlc[ix])
                continue
            elif len(word) >= 2 and word[0] == "ATOMIC_POSITIONS":
                Fat = True
                continue

            if Fat == True:
                try:
                    ap.append([float(word[1]), float(word[2]), float(word[3])])
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

        self.lc = self.lc * factor
        self.ap = np.array(ap)

        self.wrap_to_cell()

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
            self.lc[ix, :] = [float(x) for x in cell.find(f"a{ix+1}").text.split()]

        atoms = tree.getroot().find("output").find("atomic_structure").find("atomic_positions").findall("atom")
        ap = []
        for atom in atoms:
            atom_name = atom.get("name")
            if len(self.atomtype) == 0 or self.atomtype[-1] != atom_name:
                self.atomtype.append(atom_name)
                self.Naint.append(1)
            else:
                self.Naint[-1] += 1
            ap.append([float(ap0) for ap0 in atom.text.split()])

        # convert Cartesian to direct
        self.ap = np.array(ap) @ np.linalg.inv(self.lc)
        self.lc = self.lc * bohr

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

        ix0 = 0
        ap = []
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
                    for ix1 in range(3):
                        self.lc[ix0, ix1] = float(word[ix1+1])
                    ix0 += 1
                elif len(word) >= 2 and word[0] == "volume":
                    volume = float(word[1])
            elif Fat == True:
                if len(word) >= 2 and word[0] == "newtype":
                    self.Ntype += 1
                    self.atomtype.append(word[1])
                    self.Naint.append(0)
                elif len(word) >= 4 and word[0] == "coord":
                    ap.append([float(word[1]), float(word[2]), float(word[3])])
                    self.Naint[-1] += 1

        self.ap = np.array(ap)
        self.Natom = len(self.ap)

        detlc = self.volume()
        factor = (volume/detlc)**(1.0/3.0)*bohr
        self.lc = self.lc * factor

        self.wrap_to_cell()

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

        reading_ap = False
        reading_lc = False
        reading_k_grid = False
        reading_k_grid_shift = False
        factor_lc = bohr
        ap = []
        for il in range(len(line)):
            word = line[il].replace(":", " ").replace("=", " ").split()
            if len(word) == 0 or word[0][0] == "#" or word[0][0] == "!":
                continue
            if len(word) >= 2 and word[0].lower() == "begin" and word[1].lower() == "atom_coord":
                reading_ap = True
                continue
            if len(word) >= 2 and word[0].lower() == "end" and word[1].lower() == "atom_coord":
                reading_ap = False
                continue
            if len(word) >= 2 and word[0].lower() == "begin" and word[1].lower() == "cell_shape":
                reading_lc = True
                ilc = 0
                continue
            if len(word) >= 2 and word[0].lower() == "end" and word[1].lower() == "cell_shape":
                reading_lc = False
                continue
            if len(word) >= 2 and word[0].lower() == "begin" and word[1].lower() == "monkhorst_pack_grid":
                reading_k_grid = True
                continue
            if len(word) >= 2 and word[0].lower() == "end" and word[1].lower() == "monkhorst_pack_grid":
                reading_k_grid = False
                continue
            if len(word) >= 2 and word[0].lower() == "begin" and word[1].lower() == "monkhorst_pack_shift":
                reading_k_grid_shift = True
                continue
            if len(word) >= 2 and word[0].lower() == "end" and word[1].lower() == "monkhorst_pack_shift":
                reading_k_grid_shift = False
                continue
            if word[0].lower() == "boundary_sphere_radius":
                alat = float(word[1]) * 2.0  # in Bohr here
                if self.Ndim == 0:
                    np.fill_diagonal(self.lc, alat)
                elif self.Ndim == 2:
                    self.lc[2, 2] = alat
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
            if reading_lc:
                for ix in range(self.Ndim):
                    self.lc[ilc, ix] = float(word[ix])
                ilc += 1
            if reading_ap:
                ap.append([float(word[0]), float(word[1]), float(word[2])])
                self.Naint[-1] += 1
                self.Natom += 1
            if reading_k_grid:
                self.k_grid = [1, 1, 1]
                for iw in range(min(len(word), 3)):
                    self.k_grid[iw] = int(word[iw])
            if reading_k_grid_shift:
                self.k_grid_shift = [0.0, 0.0, 0.0]
                for iw in range(min(len(word), 3)):
                    self.k_grid_shift[iw] = float(word[iw])

        self.lc = self.lc * factor_lc
        self.ap = np.array(ap) * factor_ap

        if lcart:
            self.ap = self.ap @ np.linalg.inv(self.lc)
        if self.Ndim == 0:
            shift = np.array([0.5, 0.5, 0.5])
        elif self.Ndim == 2:
            shift = np.array([0.0, 0.0, 0.5])
        else:  # bulk
            shift = np.array([0.0, 0.0, 0.0])
        self.ap += shift

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
        ap = []
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
            ap.append([float(word[1]), float(word[2]), float(word[3])])

        self.ap = np.array(ap)

        if lbohr:
            self.ap = self.ap * bohr

        radius = np.max(np.linalg.norm(self.ap, axis=1))
        lc0 = 2.0 * (radius+5.0)  # vacuum is set to 5 A
        np.fill_diagonal(self.lc, lc0)

        self.ap = self.ap / lc0 + 0.5

########################################################################

    def write_vasp(self, filename="POSCAR.new"):
        f1 = open(filename, "w")
        f1.write(self.title+'\n')
        f1.write('1.0\n')
        for ix0 in range(3):
            for ix1 in range(3):
                f1.write(f"{self.lc[ix0, ix1]:21.12f}")
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
        for ia in range(self.Natom):
            for ix in range(3):
                f1.write(f"{self.ap[ia, ix]:18.12f}")
            if self.seldyn is not None:
                for ix in range(3):
                    if self.seldyn[ia][ix]:
                        f1.write(" T")
                    else:
                        f1.write(" F")
            f1.write('\n')
        f1.close()

    def write_qe(self, filename="qe_st.dat"):
        mass = load_atom_mass()

        f2 = open(filename, "w")
        f2.write("CELL_PARAMETERS angstrom\n")
        for ix0 in range(3):
            for ix1 in range(3):
                f2.write(f"{self.lc[ix0, ix1]:18.12f}")
            f2.write("\n")

        f2.write("ATOMIC_SPECIES\n")
        for i in range(self.Ntype):
            f2.write("  " + self.atomtype[i] + "  " + str(mass[self.atomtype[i]]) + "  " + self.atomtype[i] + ".upf\n")

        f2.write("ATOMIC_POSITIONS crystal\n")
        itype = 0
        iaint = 0
        for ia in range(self.Natom):
            f2.write(f"  {self.atomtype[itype]:2s}")
            for ix in range(3):
                f2.write(f"{self.ap[ia, ix]:18.12f}")
            f2.write("\n")
            iaint = iaint + 1
            if iaint == self.Naint[itype]:
                itype = itype + 1
                iaint = 0

        k_grid = self.find_k_grid()
        f2.write("K_POINTS automatic\n")
        f2.write(f"  {k_grid[0]:d}  {k_grid[1]:d}  {k_grid[2]:d} 0 0 0\n")
        f2.close()

    def write_prt(self, filename="prt_st.dat"):
        bohr = load_constant("bohr")

        volume = self.volume()*bohr**(-3)
        # convert Angstrom^3 to bohr^3

        f2 = open(filename, "w")
        f2.write("begin latticevecs\n")
        for ix0 in range(3):
            f2.write("coord")
            for ix1 in range(3):
                f2.write(f"{self.lc[ix0, ix1]:18.12f}")
            f2.write("\n")
        f2.write(f"volume {volume:24.16f}\nend latticevecs\n\nbegin coordinates\n")
        ia = 0
        for itype in range(self.Ntype):
            f2.write("newtype "+self.atomtype[itype]+"\n")
            for _ in range(self.Naint[itype]):
                f2.write("coord")
                for ix in range(3):
                    f2.write(f"{self.ap[ia, ix]:18.12f}")
                f2.write("\n")
                ia += 1
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
            radius = max(self.lc[0, 0], self.lc[1, 1], self.lc[2, 2])*0.5/funit
            if lbohr:
                f2.write(f"boundary_sphere_radius {radius:.12g}\n\n")
            else:
                f2.write(f"boundary_sphere_radius {radius:.12g} ang\n\n")
        elif self.Ndim == 1:
            f2.write("boundary_conditions wire\n")
            radius = max(self.lc[1, 1], self.lc[2, 2])*0.5/funit
            if lbohr:
                f2.write(f"boundary_sphere_radius {radius:.12g}\n\n")
            else:
                f2.write(f"boundary_sphere_radius {radius:.12g} ang\n\n")
        elif self.Ndim == 2:
            f2.write("boundary_conditions slab\n")
            radius = self.lc[2, 2]*0.5/funit
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
            for ix0 in range(self.Ndim):
                for ix1 in range(3):
                    f2.write(f"{self.lc[ix0, ix1]/funit:18.12f}")
                f2.write("\n")
            f2.write("end cell_shape\n\n")

            f2.write("kpoint_method mp\n\n")
            f2.write("begin monkhorst_pack_grid\n")
            k_grid = self.find_k_grid()
            for ix in range(self.Ndim):
                f2.write(f"  {k_grid[ix]:d}")
            f2.write("\nend monkhorst_pack_grid\n\n")
            k_grid_shift = self.find_k_grid_shift()
            f2.write("begin monkhorst_pack_shift\n")
            for ix in range(self.Ndim):
                f2.write(f"  {k_grid_shift[ix]:f}")
            f2.write("\nend monkhorst_pack_shift\n\n")

        f2.write("atom_types_num " + str(len(self.atomtype)) + "\n")
        if lcartesian:
            if lbohr:
                f2.write("coordinate_unit cartesian_bohr\n\n")
            else:
                f2.write("coordinate_unit cartesian_ang\n\n")
            apc = self.cartesian(factor=1.0/funit)
        else:
            f2.write("coordinate_unit lattice_vectors\n\n")
        ia = 0
        for itype in range(self.Ntype):
            f2.write("atom_type " + self.atomtype[itype] + "\n")
            # f2.write("Pseudopotential_Format: \n")
            # f2.write("Core_Cutoff_Radius: \n")
            f2.write("local_component s\n")  # to be modified
            # f2.write("Potential_Num: \n\n")
            # f2.write("begin Electron_Per_Orbital\n")
            # f2.write("# S P D F\n \n")
            # f2.write("end Electron_Per_Orbital\n\n")

            f2.write("begin atom_coord\n")
            if lcartesian:
                for _ in range(self.Naint[itype]):
                    for ix in range(3):
                        f2.write(f"{apc[ia, ix]:18.12f}")
                    f2.write("\n")
                    ia += 1
            else:
                for _ in range(self.Naint[itype]):
                    for ix in range(3):
                        f2.write(f"{self.ap[ia, ix]:18.12f}")
                    f2.write("\n")
                    ia += 1
            f2.write("end atom_coord\n\n")

        f2.close()

    def write_wannier90(self, filename="wannier90_st.dat"):
        f2 = open(filename, "w")
        f2.write("Begin Unit_Cell_Cart\n")
        for ix0 in range(3):
            for ix1 in range(3):
                f2.write(f"{self.lc[ix0, ix1]:18.12f}")
            f2.write("\n")
        f2.write("End Unit_Cell_Cart\n\n")

        f2.write("Begin Projections\n")
        f2.write("  random\n")
        for i in range(self.Ntype):
            f2.write("  "+self.atomtype[i]+"  :\n")
        f2.write("End Projections\n\n")

        f2.write("Begin Atoms_Frac\n")
        ia = 0
        for itype in range(self.Ntype):
            for _ in range(self.Naint[itype]):
                f2.write(f"  {self.atomtype[itype]:2s}")
                for ix in range(3):
                    f2.write("{self.ap[ia, ix]:18.12f}")
                f2.write("\n")
                ia += 1
        f2.write("End Atoms_Frac\n\n")
        f2.close()

    def write_xyz(self, filename="structure.xyz"):

        self.Ndim = 0

        f2 = open(filename, "w")
        f2.write(str(self.Natom) + "\n\n")
        ia = 0
        apc = self.cartesian()
        for itype in range(self.Ntype):
            for _ in range(self.Naint[itype]):
                f2.write(f"  {self.atomtype[itype]:2s}")
                for ix in range(3):
                    f2.write(f"{apc[ia, ix]:18.12f}")
                f2.write("\n")
                ia += 1

        f2.close()

########################################################################

    def wrap_to_cell(self):
        tol = 1e-8
        self.ap = self.ap - self.ap//1.0
        self.ap[(self.ap < tol) | (self.ap > 1.0-tol)] = 0.0

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
            disp = np.array(disp) @ rlc
        self.ap = self.ap + disp
        if lbox:
            self.wrap_to_cell()

    def rotate(self, theta):
        # rotate the structure around the z-axis (in degree)
        t = math.radians(theta)
        c = math.cos(t)
        s = math.sin(t)
        M = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
        self.lc = self.lc @ M

    def flip(self):
        # flip the z coordinate of the structure
        self.ap[:, 2] = 1.0 - self.ap[:, 2]
        self.wrap_to_cell()

    def supercell(self, N):
        # N should be a list of 3 integers
        N_total = N[0]*N[1]*N[2]

        M = np.array([[N[0], 0, 0], [0, N[1], 0], [0, 0, N[2]]])
        self.lc = M @ self.lc

        iz, iy, ix = np.meshgrid(np.arange(N[2]), np.arange(N[1]), np.arange(N[0]), indexing='ij')
        shift = np.stack([ix, iy, iz], axis=-1).reshape(-1, 3) / np.array(N)

        ap = np.zeros((N_total*self.Natom, 3))
        for ia in range(self.Natom):
            ap[ia*N_total:(ia+1)*N_total, :] = self.ap[ia, :]/np.array(N)+shift
        self.ap = ap

        for itype in range(self.Ntype):
            self.Naint[itype] = self.Naint[itype]*N_total
        self.Natom = self.Natom*N_total

        self.wrap_to_cell()

    def vacuum(self, z_vac):
        # Add a vacuum layer to the structure
        # z_vac is in angstrom
        # This version assumes a1 x a2 // a3 // z
        if z_vac <= 0:
            print("Error: z_vac <= 0")
            sys.exit()

        a3_new = self.lc[2, 2] + z_vac

        self.wrap_to_cell()

        ia = 0
        ap_new = []
        for itype in range(self.Ntype):
            ap_new0 = []
            for iaint in range(self.Naint[itype]):
                if self.ap[ia, 2] < 1e-8:
                    ap_new0.append([self.ap[ia, 0], self.ap[ia, 1], 1.0])
                ia += 1
            ap_new.append(ap_new0)

        for itype in range(self.Ntype):
            for ap in ap_new[itype]:
                self.add_atom(itype, ap)

        self.ap[:, 2] = (self.ap[:, 2] - 0.5) * self.lc[2, 2] / a3_new + 0.5

        self.lc[2, 2] = a3_new

    def add_atom(self, itype, ap, new_type=None, add_to_head=False):
        # If new_type is defined, a new atom type will be inserted.
        # add_to_head defines if the atom will be added to the head.
        if itype < 0 or itype > self.Ntype:
            print("Error: itype ("+str(itype)+") is not in range [0, "+str(self.Ntype)+"] in add_atom function.")
            sys.exit()

        if new_type is None and itype == self.Ntype:
            print("Warning: itype ("+str(itype)+") == self.Ntype in add_atom function, but atom name is not defined.")
            new_type = "X"

        if new_type is not None:  # This is a new type
            self.Ntype += 1
            self.atomtype.insert(itype, new_type)
            self.Naint.insert(itype, 0)

        self.Natom += 1
        if add_to_head:
            insert_index = sum(self.Naint[:itype])
        else:
            insert_index = sum(self.Naint[:itype+1])
        self.ap = np.insert(self.ap, insert_index, ap, axis=0)
        self.Naint[itype] += 1

    def delete_atom(self, ia):
        if ia < 0 or ia >= self.Natom:
            print("Error: ia ("+ia+") is out of range in delete_atom function")
            sys.exit()

        self.ap = np.delete(self.ap, ia, axis=0)
        self.Natom -= 1

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

        ap = self.ap[ia, :]
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
        apc = self.cartesian(shift=0.0)
        self.ap = apc @ np.linalg.inv(lc_new)
        self.lc = lc_new

    def rlc(self):
        pi = load_constant("pi")
        return np.linalg.inv(self.lc) * (2.0*pi)

    def volume(self):
        return abs(np.linalg.det(self.lc))

    def cartesian(self, shift=None, factor=1.0):
        if shift is None:
            if self.Ndim == 3:
                shift = np.array([0.0, 0.0, 0.0])
            elif self.Ndim == 2:
                shift = np.array([0.0, 0.0, -0.5])
            elif self.Ndim == 0:
                shift = np.array([-0.5, -0.5, -0.5])
        return (self.ap + shift) @ self.lc * factor

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
                        atom_flag[int(atom_list[0]) - 1] = True
                    else:
                        atom_flag[int(atom_list[0]) - 1: int(atom_list[1])] = True
        return atom_flag

    def find_k_grid(self):
        if self.k_grid is not None:
            return self.k_grid
        grid = [1] * 3
        for ix in range(self.Ndim):
            grid[ix] = math.ceil(30.0/(np.linalg.norm(self.lc[ix])))
        return grid

    def find_k_grid_shift(self):
        if self.k_grid_shift is not None:
            return self.k_grid_shift
        return [0.0] * 3

    def find_ndim(self):
        # if c not vertical: 3D
        # if atoms span nearly whole c: 3D
        # else if a has nonzero y: 2D
        # else if atoms span nearly whole b: 2D
        # else if atoms span nearly whole a: 1D
        # else: 0D

        if self.Ndim != 3:
            print("Warning: find_ndim called, but Ndim != 3")
            return
        vacuum = 5.0

        if not np.allclose(self.lc[:2, 2], 0.0):
            self.Ndim = 3
            return

        apc = self.cartesian()
        if max(apc[:, 2]) > self.lc[2, 2] - vacuum or min(apc[:, 2]) < vacuum:
            self.Ndim = 3
            return

        if not np.isclose(self.lc[0, 1], 0.0):
            self.Ndim = 2

        if max(apc[:, 1]) > self.lc[1, 1] - vacuum or min(apc[:, 1]) < vacuum:
            self.Ndim = 2
            return

        if max(apc[:, 0]) > self.lc[0, 0] - vacuum or min(apc[:, 0]) < vacuum:
            self.Ndim = 1
            return

        self.Ndim = 0
        return
