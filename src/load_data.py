import os
def load_package_name() :
    # read the package name from the data/package_name.dat into a dictionary
    # data/package_name.dat: 1st column is the key, following columns are the values
    this_dir, this_file = os.path.split(__file__)
    data_path = os.path.join(this_dir, "..", "data", "package_name.dat")
    f0=open(data_path,"r")
    line=f0.readlines()
    f0.close()
    package_name={}
    for l in line :
        word=l.split()
        if len(word)==0 or word[0][0] in {"#", "!"} :
            continue
        package_name[word[0]]=word[1:]
    return package_name

def load_palette() :
    # read the color palette from the data/palette.dat into a dictionary :
    # data/palette.dat: key, values
    this_dir, this_file = os.path.split(__file__)
    data_path = os.path.join(this_dir, "..", "data", "palette.dat")
    f0=open(data_path,"r")
    line=f0.readlines()
    f0.close()
    palette={}
    for l in line :
        word=l.split()
        if len(word)==0 or word[0][0] in {"#", "!"} :
            continue
        palette[word[0]]="#"+word[1]
    return palette

def load_symops() :
    # read the crystallographic symmetry operations from the data/symops.dat into a list of length 48 :
    # data/symops.dat: index, xyz, opname, mtx[9]
    this_dir, this_file = os.path.split(__file__)
    data_path = os.path.join(this_dir, "..", "data", "symops.dat")
    f0=open(data_path,"r")
    line=f0.readlines()
    f0.close()
    mtx=[]
    name=[]
    for l in line :
        word=l.split()
        if len(word)==0 or word[0][0] in {"#", "!"} :
            continue
        mtx0=[]
        for ix in range(3):
            mtx0.append([int(word[ix*3+3]),int(word[ix*3+4]),int(word[ix*3+5])])
        mtx.append(mtx0)
        name.append(word[2])
    return mtx, name

def load_constant(const_name) :
    # read the constants from the data/constant.dat file
    this_dir, this_file = os.path.split(__file__)
    data_path = os.path.join(this_dir, "..", "data", "constant.dat")
    f0=open(data_path,"r")
    line=f0.readlines()
    f0.close()
    const_name=const_name.lower()
    for l in line :
        word=l.split()
        if len(word)==0 or word[0][0] in {"#", "!"} :
            continue
        if word[0]==const_name :
            const=float(word[1])
    return const

def load_atom_index() :
    # return a dictionary of atom_name -> atom_index
    this_dir, this_file = os.path.split(__file__)
    data_path = os.path.join(this_dir, "..", "data", "atom.dat")
    f0=open(data_path,"r")
    line=f0.readlines()
    f0.close()

    index={}
    for l in line :
        word=l.split()
        if len(word)==0 or word[0][0] in {"#", "!"} :
            continue
        index[word[2]]=int(word[0])
    return index

def load_atom_mass() :
    # return a dictionary of atom_name -> atom_mass
    this_dir, this_file = os.path.split(__file__)
    data_path = os.path.join(this_dir, "..", "data", "atom.dat")
    f0=open(data_path,"r")
    line=f0.readlines()
    f0.close()

    mass={}
    for l in line:
        word=l.split()
        if len(word)==0 or word[0][0] in {"#", "!"} :
            continue
        mass[word[2]]=float(word[3])
    return mass

def load_atom_color():
    # return a dictionary of atom_name -> atom_color
    this_dir, this_file = os.path.split(__file__)
    data_path = os.path.join(this_dir, "..", "data", "atom.dat")
    f0=open(data_path,"r")
    line=f0.readlines()
    f0.close()

    color={}
    for l in line:
        word=l.split()
        if len(word)==0 or word[0][0] in {"#", "!"} :
            continue
        color[word[2]]=word[4]
    return color

