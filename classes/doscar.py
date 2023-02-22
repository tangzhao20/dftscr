import xml.etree.ElementTree as ET

class DOSCAR :
    # ef : float
    # ispin : integer
    # nedos : integer
    def __init__(self, filename="DOSCAR",empty=False) :
        if empty:
            return
        f0=open(filename,"r")
        line=f0.readlines()
        f0.close()
        self.ef=float(line[5].split()[3])
        del line

    def fileread_xml(self, filename="pwscf.xml") :
        tree=ET.parse(filename)
        self.ef=float(tree.getroot().find("output").find("band_structure").find("fermi_energy").text)*27.211386245988

