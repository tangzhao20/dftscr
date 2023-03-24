import sys
import xml.etree.ElementTree as ET

class DOSCAR :
    # ef : float
    # Ns : integer
    # Nedos : integer
    # energy[Nedos]
    # dos[Nedos]
    def __init__(self, filename="DOSCAR",empty=False) :
        if empty:
            return
        f0=open(filename,"r")
        line=f0.readlines()
        f0.close()
        if len(line)<7 :
            print("Error: incomplete DOSCAR")
            sys.exit()
        self.Nedos=int(line[5].split()[2])
        self.ef=float(line[5].split()[3])
        self.Ns=(len(line[6].split())-1)//2
        self.energy=[0.0]*self.Nedos
        if self.Ns==1 :
            self.dos=[[0.0]*self.Nedos]
        else :
            self.dos=[[0.0]*self.Nedos,[0.0]*self.Nedos]

        for il in range(self.Nedos) :
            word=line[il+6].split()
            self.energy[il]=float(word[0])
            self.dos[0][il]=float(word[1])
            if self.Ns==2 :
                self.dos[1][il]=float(word[2])
            
        del line

    def fileread_xml(self, filename="pwscf.xml") :
        tree=ET.parse(filename)
        self.ef=float(tree.getroot().find("output").find("band_structure").find("fermi_energy").text)*27.211386245988


    def energyshift(self,ezero) :
        for ie in range(self.Nedos) :
            self.energy[ie]=self.energy[ie]-ezero
