import numpy as np
from load_data import load_atom_mass


def match(poscar1, poscar2):
    if np.linalg.norm(poscar1.lc-poscar2.lc) > 1e-6:
        print("Lattice vectors don't match.")
        return False
    if poscar1.Natom != poscar2.Natom or poscar1.Ntype != poscar2.Ntype:
        print("Atoms don't match: #atoms or #types")
        return False
    for i in range(poscar1.Ntype):
        if poscar1.atomtype[i] != poscar2.atomtype[i] or poscar1.Naint[i] != poscar2.Naint[i]:
            print("Atoms don't match: "+str(i+1)+": " +
                  poscar1.atomtype[i]+" "+str(poscar1.Naint[i])+" vs "+poscar2.atomtype[i]+" "+str(poscar2.Naint[i]))
            return False
    print("Structures match each other.")
    return True


def moveatoms(poscar1, poscar2, factor):
    # move atoms of poscar1 towards a new structure poscar2
    for i in range(poscar1.Natom):
        for j in range(3):
            if poscar2.ap[i][j]-poscar1.ap[i][j] >= 0.5:
                newap = poscar2.ap[i][j]-1.0
            elif poscar2.ap[i][j]-poscar1.ap[i][j] < -0.5:
                newap = poscar2.ap[i][j]+1.0
            else:
                newap = poscar2.ap[i][j]
            poscar1.ap[i][j] = poscar1.ap[i][j]*(1.0-factor)+newap*factor
    poscar1.wrap_to_cell()


def displacement(poscar1, poscar2):
    # find the displacement vector of two POSCARs
    # disp[Na][3]
    disp = np.array(poscar2.ap)-np.array(poscar1.ap)
    disp = (disp+0.5) % 1-0.5
    disp = disp @ poscar1.lc
    disp = disp.tolist()
    return disp


def general_q(poscar1, poscar2):
    # calculate the generalized coordinate of poscar2 compared to poscar1

    disp = displacement(poscar1, poscar2)

    mass = load_atom_mass()

    mass_list = []
    for it in range(poscar1.Ntype):
        for ia in range(poscar1.Naint[it]):
            mass_list.append(mass[poscar1.atomtype[it]])
    mass_list = np.sqrt(mass_list)

    Q = np.dot(np.linalg.norm(disp, axis=1),  mass_list)
    return Q
