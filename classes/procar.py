import sys
import numpy as np
import re

class PROCAR :
    # Nk
    # Nb
    # Na
    # Nlm
    # proj[k][b][a][orb]
    # orbname[orb]
 
    def __init__(self,filename="PROCAR",empty=False) :
        if empty :
            self.Nk=0
            self.Nb=0
            self.Na=0
            self.Nlm=0
            self.proj=[]
            self.orb=[]
            self.orbname=[]
            return
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
        a=-1
        for il in range(len(line)) :
            word=line[il].split()
            if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
                continue
            if word[0]=="k-point" :
                k=int(word[1])
                proj0=[]
                #word=line[].split()
            elif word[0]=="band" :
                proj1=[]
                b=int(word[1])
            elif word[0]=="ion" :
                a=0
            elif a>=0 and a<self.Na :
                if int(word[0])!=a+1 :
                    print("atom number mismatch on line "+str(il+1))
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

        if self.Nlm==9 :
            self.orbname=["s","py","pz","px","dxy","dyz","dz2","dxz","x2-y2"]
        elif self.Nlm==16 :
            self.orbname=["s","py","pz","px","dxy","dyz","dz2","dxz","x2-y2","fy3x2","fxyz","fyz2","fz3","fxz2","fzx2","fx3"]
        else :
            print("Nlm="+str(self.Nlm)+" is not support yet")
            sys.exit()

        del proj0
        del proj1
        del proj2
        del line
        del word

    def fileread_qe(self, filename="projwfc.out") :
        f0=open(filename,"r")
        line=f0.readlines()
        f0.close()

        lmmap=[]
        atommap=[]
        fread=False
        ffirstk=True
        for il in range(len(line)) :
            word=line[il].split()
            if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
                continue
            if word[0]=="state" :
                ia=int(word[4])-1
                if len(atommap)==0 or ia!=atommap[-1] :
                    ilm=0
                else :
                    ilm+=1
                lmmap.append(ilm)
                atommap.append(ia)
                self.Nlm=max(ilm+1,self.Nlm)
                self.Na=max(ia+1,self.Na)
            elif word[0]=="nkstot" :
                self.Nk=int(word[2])
            elif word[0]=="nbnd" :
                self.Nb=int(word[2])
            elif word[0]=="k" :
                if ffirstk:
                    ffirstk=False
                    for ik in range(self.Nk) :
                        proj0=[]
                        for ib in range(self.Nb) :
                            proj1=[]
                            for ia in range(self.Na) :
                                proj1.append([0.0]*self.Nlm)
                            proj0.append(proj1)
                            del proj1
                        self.proj.append(proj0)
                        del proj0
                    ik=-1
                ik+=1
                ib=-1
            elif word[0]=="psi" :
                ib+=1
                fread=True
            elif word[0]=="|psi|^2" :
                fread=False

            if fread :
                number=re.findall(r"[-+]?(?:\d*\.*\d+)", line[il]) # find all numbers
                for ii in range(0,len(number),2) :
                    # proj[k][b][a][orb]
                    self.proj[ik][ib][atommap[int(number[ii+1])-1]][lmmap[int(number[ii+1])-1]]=float(number[ii])
        if self.Nlm==4 :
            self.orbname=["s","pz","px","py"]
        elif self.Nlm==9 :
            self.orbname=["s","pz","px","py","dz2","dxz","dyz","x2-y2","dxy"]
        elif self.Nlm==16 :
            self.orbname=["s","pz","px","py","dz2","dxz","dyz","x2-y2","dxy","fz3","fxz2","fyz2","fzx2","fxyz","fx3","fy3x2"]
        else :
            print("Nlm="+str(self.Nlm)+" is not support yet")
            sys.exit()
                
    def plot(self,atomflag,orbflag,fac) :
        plotproj=[]
        for b in range(self.Nb) :
            plotproj0=[]
            for k in range(self.Nk) :
                plotproj1=0.0
                for a in range(self.Na) :
                    for i in range(self.Nlm) :
                        if atomflag[a] and orbflag[i]==1 :
                            plotproj1+=self.proj[k][b][a][i]
                plotproj0.append(plotproj1*fac)
            plotproj.append(plotproj0)
        return plotproj

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
