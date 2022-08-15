import numpy as np

class PROCAR :
    def __init__(self,filename="PROCAR") :
        f0=open(filename,"r")
        line=f0.readlines()
        f0.close()
        word=line[1].split()
        self.Nk=int(word[3])
        self.Nb=int(word[7])
        self.Na=int(word[11])

        word=line[7].split()
        self.Nlm=len(word)-2

        print("Na "+str(self.Na)+" Nk "+str(self.Nk)+" Nb "+str(self.Nb)+" Nlm "+str(self.Nlm))
            
        self.proj=[]
# proj[k][b][a][orb]
# orb: s py pz px dxy dyz dz2 dxz dx2-y2
        a=-1
        for l in range(len(line)) :
            word=line[l].split()
            if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
                continue
            if word[0]=="k-point" :
                k=int(word[1])
                proj0=[]
    #            word=line[].split()
            elif word[0]=="band" :
                proj1=[]
                b=int(word[1])
            elif word[0]=="ion" :
                a=0
            elif a>=0 and a<self.Na :
                if int(word[0])!=a+1 :
                    print("atom number mismatch on line "+str(l+1))
                    sys.exit()
                proj2=[]
                for lm in range(self.Nlm) :
                    proj2.append(float(word[lm+1]))
                proj1.append(proj2)
                if a<self.Na-1 :
                    a=a+1
                else :
                    a=-1
                    proj0.append(proj1)
                    if b==self.Nb-1 :
                        self.proj.append(proj0)
        del proj0
        del proj1
        del proj2
        del line
        del word

    def plot(self,atomflag,orbflag,fac) :
        plotproj=[]
        for b in range(self.Nb) :
            plotproj0=[]
            for k in range(self.Nk) :
                plotproj1=0.0
                for a in range(self.Na) :
                    for i in range(self.Nlm) :
                        if atomflag[a] and orbflag[i]==1 :
                            plotproj1=plotproj1+self.proj[k][b][a][i]
                plotproj0.append(plotproj1*fac)
            plotproj.append(plotproj0)
        return plotproj

    def readorblist(self, orblist) :
        if self.Nlm==9 :
            orbname=["s","py","pz","px","dxy","dyz","dz2","dxz","x2-y2"]
        elif self.Nlm==16 :
            orbname=["s","py","pz","px","dxy","dyz","dz2","dxz","x2-y2","fy3x2","fxyz","fyz2","fz3","fxz2","fzx2","fx3"]
        orblist0=orblist.split("+")
        orbflag=[]
        for i in range(self.Nlm) :
            orbflag.append(0)
        for i in range(len(orblist0)) :
            Fprojvalid=False
            for j in range(self.Nlm) :
                if orblist0[i]==orbname[j] :
                    orbflag[j]=1
                    Fprojvalid=True
                    break
            if Fprojvalid==True :
                continue
            if orblist0[i]=="p" :
                for j in range(1,4):
                    orbflag[j]=1
                Fprojvalid=True
            elif orblist0[i]=="d" :
                for j in range(4,9):
                    orbflag[j]=1
                Fprojvalid=True
            elif orblist0[i]=="f" :
                if self.Nlm==16 :
                    for j in range(9,16):
                        orbflag[j]=1
                    Fprojvalid=True
            elif orblist0[i]=="dx2-y2" :
                orbflag[8]=1
                Fprojvalid=True
            elif orblist0[i]=="all" :
                for j in range(self.Nlm):
                    orbflag[j]=1
                Fprojvalid=True
            if Fprojvalid==False :
                print("projector "+orblist0[i]+" does not exist")
        return orbflag

    def readatomlist(self, atomlist, poscar1) :
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
                for j in range(len(poscar1.atomtype)) :
                    Iatom1=Iatom1+poscar1.Naint[j]
                    if atomlist0[i]==poscar1.atomtype[j] :
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
