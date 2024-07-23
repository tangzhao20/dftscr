import sys
import numpy as np
import re

class Procar:
    # Nk
    # Nb
    # Na
    # Nlm
    # proj[k][b][a][orb]
    # orb_name[orb]
 
    def __init__(self):
        self.Nk=0
        self.Nb=0
        self.Na=0
        self.Nlm=0
        self.proj=[]
        self.orb=[]
        self.orb_name=[]

    def read_vasp(self, filename="PROCAR"):

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
        for il in range(len(line)):
            word=line[il].split()
            if len(word)==0 or word[0][0]=="#" or word[0][0]=="!":
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
            self.orb_name=["s","py","pz","px","dxy","dyz","dz2","dxz","x2-y2"]
        elif self.Nlm==16 :
            self.orb_name=["s","py","pz","px","dxy","dyz","dz2","dxz","x2-y2","fy3x2","fxyz","fyz2","fz3","fxz2","fzx2","fx3"]
        else :
            print("Nlm="+str(self.Nlm)+" is not support yet")
            sys.exit()

        del proj0
        del proj1
        del proj2
        del line
        del word

    def read_qe(self, filename="projwfc.out"):
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
            self.orb_name=["s","pz","px","py"]
        elif self.Nlm==9 :
            self.orb_name=["s","pz","px","py","dz2","dxz","dyz","x2-y2","dxy"]
        elif self.Nlm==16 :
            self.orb_name=["s","pz","px","py","dz2","dxz","dyz","x2-y2","dxy","fz3","fxz2","fyz2","fzx2","fxyz","fx3","fy3x2"]
        else :
            print("Nlm="+str(self.Nlm)+" is not support yet")
            sys.exit()
                
    def plot(self,atom_flag,orb_flag,fac):
        plot_proj=[]
        for b in range(self.Nb) :
            plot_proj0=[]
            for k in range(self.Nk) :
                plot_proj1=0.0
                for a in range(self.Na) :
                    for i in range(self.Nlm) :
                        if atom_flag[a] and orb_flag[i]==1 :
                            plot_proj1+=self.proj[k][b][a][i]
                plot_proj0.append(plot_proj1*fac)
            plot_proj.append(plot_proj0)
        return plot_proj

    def read_orb_list(self, orb_list):
        orb_list0=orb_list.split("+")
        orb_list1=[]
        orb_flag=[0]*self.Nlm
        for j in orb_list0 :
            if j in self.orb_name:
                orb_list1.append(j)
            elif j=="p" :
                orb_list1+=["px","py","pz"]
            elif j=="d" :
                orb_list1+=["dxy","dyz","dz2","dxz","x2-y2"]
            elif j=="f" :
                orb_list1+=["fy3x2","fxyz","fyz2","fz3","fxz2","fzx2","fx3"]
            elif j=="dx2-y2" :
                orb_list1.append("x2-y2")
            elif j=="all" :
                orb_list1+=self.orb_name
            else :
                print("projector "+j+" does not exist")

        for j in orb_list1:
            if j in self.orb_name:
                orb_flag[self.orb_name.index(j)]=1
            else :
                print("projector "+j+" does not exist")

        return orb_flag

    def read_atom_list(self, atom_list, poscar1):
        atom_flag=[]
        if atom_list=="all" :
            for a in range(self.Na) :
                atom_flag.append(1)
        else :
            for a in range(self.Na) :
                atom_flag.append(0)
            atom_list0=atom_list.split(",")
            for i in range(len(atom_list0)) :
                Fatom_valid=False
                Iatom0=0
                Iatom1=0
                for j in range(len(poscar1.atomtype)) :
                    Iatom1=Iatom1+poscar1.Naint[j]
                    if atom_list0[i]==poscar1.atomtype[j] :
                        for a in range(Iatom0,Iatom1) :
                            atom_flag[a]=1
                        Fatom_valid=True
                    Iatom0=Iatom1
                if Fatom_valid==False :
                    atom_list1=atom_list0[i].split("-")
                    if len(atom_list1)==1 :
                        atom_flag[int(atom_list1[0])-1]=1
                        Fatom_valid=True
                    else :
                        for j in range(int(atom_list1[0])-1,int(atom_list1[1])) :
                            atom_flag[j]=1
                        Fatom_valid=True
        return atom_flag
