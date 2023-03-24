import os
def load_packagename():
    # read the package name from the data/packagename.dat into a dictionary
    # data/packagename.dat: 1st column is the key, following columns are the values
    this_dir, this_filename = os.path.split(__file__)
    DATA_PATH = os.path.join(this_dir, "..", "data", "packagename.dat")
    f0=open(DATA_PATH,"r")
    line=f0.readlines()
    f0.close()
    packagename={}
    for l in line :
        word=l.split()
        if len(word)==0 or word[0][0] in {"#", "!"} :
            continue
        packagename[word[0]]=word[1:]
    return packagename

def load_palette() :
    # read the color palette from the data/palette.dat into a list [blue, orange, gray, white, black] :
    # data/palette.dat: key, values
    this_dir, this_filename = os.path.split(__file__)
    DATA_PATH = os.path.join(this_dir, "..", "data", "palette.dat")
    f0=open(DATA_PATH,"r")
    line=f0.readlines()
    f0.close()
    palette=[]
    for l in line :
        word=l.split()
        if len(word)==0 or word[0][0] in {"#", "!"} :
            continue
        palette.append("#"+word[1])
    return palette
