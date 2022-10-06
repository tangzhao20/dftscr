class DOSCAR :
    # ef : float
    # ispin : integer
    # nedos : integer
    def __init__(self, filename="DOSCAR") :
        f0=open(filename,"r")
        line=f0.readlines()
        f0.close()
        self.ef=float(line[5].split()[3])
        del line

    def ef_out(self) :
        return(self.ef)
