import sys
import os
import re
import numpy as np
import xml.etree.ElementTree as ET

class DOSCAR :
    # ef : float
    # Ns : integer
    # Nedos : integer
    # energy[Nedos]
    # dos[Ns][Nedos]
    def __init__(self, filename="DOSCAR",empty=False) :
        if empty:
            self.ef=0.0
            self.Ns=0
            self.Nedos=0
            self.energy=[]
            self.dos=[[]]
            self.Lpdos=False
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
        #TODO: test the metal case
        if tree.getroot().find("output").find("band_structure").find("highestOccupiedLevel") is None :
            self.ef=float(tree.getroot().find("output").find("band_structure").find("fermi_energy").text)*27.211386245988
        else :
            self.ef=float(tree.getroot().find("output").find("band_structure").find("highestOccupiedLevel").text)*27.211386245988

    def fileread_qe(self, filename="pwscf.dos") :
        f0=open(filename,"r")
        line=f0.readlines()
        f0.close()

        if len(line)<2 :
            print("Error: incomplete "+filename)
            sys.exit()
        self.ef=float(line[0].split()[-2])
        self.Ns=(len(line[1].split())-1)//2
        if self.Ns==2 :
            self.dos.append([])
        for l in line :
            word=l.split()
            if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
                continue
            self.Nedos+=1
            self.energy.append(float(word[0]))
            self.dos[0].append(float(word[1]))
            if self.Ns==2 :
                self.dos[1].append(float(word[2]))

    def readpdos_qe(self) :
        # pdos[Ns][Nepdos][Na][Nlm]
        # x_pdos[Nepdos]
        # atomtype[Ntype]
        # Naint[Ntype]
        # orbname[Nlm]

        if self.Ns==2 :
            print("QE pdos has not been supported for spin polarized system yet.")
            sys.exit()

        files=os.listdir()
        pdoslist = [s for s in files if ".pdos_atm#" in s]
        filelist=[]
        for f in files:
            fileinfo = re.search(r"(.+?)\.pdos_atm#(\d+)\((.+?)\)_wfc#(\d+)\((.+?)\)", f)
            if fileinfo :
                filename=fileinfo.group(0)
                system=fileinfo.group(1)
                iatom=fileinfo.group(2)
                atomtype=fileinfo.group(3)
                iorbit=fileinfo.group(4)
                orbittype=fileinfo.group(5)
                filelist.append([filename, int(iatom), atomtype, int(iorbit)])
        if len(filelist)==0 :
            print("Error: No valid pdos file found.")
            sys.exit()
        filelist = sorted(filelist, key=lambda x: (x[1], x[3]))
        self.Na=max(filelist, key=lambda x: x[1])[1]
        self.Nlm=max(filelist, key=lambda x: x[3])[3]**2

        # Read energies and Nepdos from the first file
        f=filelist[0]
        f0=open(f[0],"r")
        line=f0.readlines()
        f0.close()
        self.energy_pdos=[]
        for l in line :
            word=l.split()
            if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
                continue
            self.energy_pdos.append(float(word[0]))
        self.Nepdos=len(self.energy_pdos)
        
        self.pdos=[]
        self.pdos=np.arange(self.Ns*self.Nepdos*self.Na*self.Nlm)
        self.pdos=np.reshape(self.pdos,(self.Ns,self.Nepdos,self.Na,self.Nlm)).tolist()

        for f in filelist :
            # f=[filename, iatom, atomtype, iorbit]
            f0=open(f[0],"r")
            line=f0.readlines()
            f0.close()
            ie=0
            for l in line :
                word=l.split()
                if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
                    continue
                if f[3]==1 :
                    self.pdos[0][ie][f[1]-1][0]=float(word[2])
                elif f[3]==2 :
                    self.pdos[0][ie][f[1]-1][1:4]=[float(word[2]),float(word[3]),float(word[4])]
                ie+=1

        # build atom list atomtype[Ntype] and Naint[Ntype]
        self.atomtype=[]
        self.Naint=[]
        for f in filelist :
            if f[3]==1 :
                if len(self.atomtype)==0 or self.atomtype[-1]!=f[2] :
                    self.atomtype.append(f[2])
                    self.Naint.append(1)
                else :
                    self.Naint[-1]+=1

        # Only s and p orbitals are supported now
        self.orbname=["s","pz","px","py"]
        self.Lpdos=True
                
################################################################################

    def energyshift(self,ezero) :
        for ie in range(self.Nedos) :
            self.energy[ie]=self.energy[ie]-ezero
        if self.Lpdos :
            for ie in range(self.Nepdos) :
                self.energy_pdos[ie]=self.energy_pdos[ie]-ezero
            
    def plotpdos(self,atomflag,orbflag) :
        pdosout=[]
        for ie in range(self.Nepdos) :
            pdosout0=0.0
            for ia in range(self.Na) :
                for im in range(self.Nlm) :
                    if atomflag[ia] and orbflag[im]==1 :
                        pdosout0+=self.pdos[0][ie][ia][im]
            pdosout.append(pdosout0)
        return pdosout

    def readorblist(self, orblist) :
        orblist0=orblist.split("+")
        orblist1=[]
        orbflag=[0]*self.Nlm
        for j in orblist0 :
            if j in self.orbname:
                orblist1.append(j)
            elif j=="p" :
                orblist1+=["px","py","pz"]
            elif j=="d" :
                orblist1+=["dxy","dyz","dz2","dxz","x2-y2"]
            elif j=="f" :
                orblist1+=["fy3x2","fxyz","fyz2","fz3","fxz2","fzx2","fx3"]
            elif j=="dx2-y2" :
                orblist1.append("x2-y2")
            elif j=="all" :
                orblist1+=self.orbname
            else :
                print("projector "+j+" does not exist")

        for j in orblist1:
            if j in self.orbname:
                orbflag[self.orbname.index(j)]=1
            else :
                print("projector "+j+" does not exist")

        return orbflag

    def readatomlist(self, atomlist) :
        atomflag=[]
        if atomlist=="all" :
            for a in range(self.Na) :
                atomflag.append(1)
        else :
            for a in range(self.Na) :
                atomflag.append(0)
            atomlist0=atomlist.split(",")
            for i in range(len(atomlist0)) :
                Fatomvalid=False
                Iatom0=0
                Iatom1=0
                for j in range(len(self.atomtype)) :
                    Iatom1=Iatom1+self.Naint[j]
                    if atomlist0[i]==self.atomtype[j] :
                        for a in range(Iatom0,Iatom1) :
                            atomflag[a]=1
                        Fatomvalid=True
                    Iatom0=Iatom1
                if Fatomvalid==False :
                    atomlist1=atomlist0[i].split("-")
                    if len(atomlist1)==1 :
                        atomflag[int(atomlist1[0])-1]=1
                        Fatomvalid=True
                    else :
                        for j in range(int(atomlist1[0])-1,int(atomlist1[1])) :
                            atomflag[j]=1
                        Fatomvalid=True
        return atomflag
