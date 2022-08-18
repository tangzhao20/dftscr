import numpy as np
import sys

class POSCAR:
    # lc[3][3]
    # Natom
    # Ntype
    # f_seldyn
    # atomtype[Ntype]
    # Naint[Ntype]
    # ap[Natom][3]
    # seldyn[Natom][3] (if f_seldyn)
    # dmass{}
    
    def __init__(self, filename="POSCAR") :
        f0=open(filename,"r")
        line=f0.readlines()
        f0.close()

        lineoff=0
        self.title=line[0]
        factor=float(line[1].split()[0])
        self.lc=[]
        for i in range(2,5) :
            word=line[i].split()
            self.lc.append([float(word[0])*factor,float(word[1])*factor,float(word[2])*factor])
        word=line[5].split()
        self.Ntype=len(word)
        self.atomtype=[]
        for i in range(self.Ntype) :
            self.atomtype.append(word[i])
        word=line[6].split()
        self.Naint=[]
        self.Natom=0
        for i in range(self.Ntype) :
            self.Naint.append(int(word[i]))
            self.Natom=self.Natom+self.Naint[i]
        word=line[7].split()
        if word[0][0]=='S' or word[0][0]=='s' :
            self.f_seldyn=True
            lineoff=lineoff+1
        else :
            self.f_seldyn=False
        word=line[7+lineoff].split()
        if word[0][0]!='D' and word[0][0]!='d' :
            print("Only atomic position 'Direct' supported")
            exit()
        self.ap=[]
        if self.f_seldyn :
            self.seldyn=[]
        for i in range(self.Natom) :
            word=line[i+lineoff+8].split()
            self.ap.append([float(word[0]),float(word[1]),float(word[2])])
            if self.f_seldyn :
                self.seldyn.append([bool(word[3]),bool(word[4]),bool(word[5])])
        self.movetobox()
        del line
        del word

        self.dmass={}
        self.dmass["B"]=10.81
        self.dmass["C"]=12.011
        self.dmass["N"]=14.007


    def movetobox(self) :
        for i in range(self.Natom) :
            for j in range(3) :
                self.ap[i][j]=self.ap[i][j]-self.ap[i][j]//1.0
                if self.ap[i][j]<0.0000000001 or self.ap[i][j]>0.9999999999 :
                    self.ap[i][j]=0.0

    def filewrite(self, filename):
        f1=open(filename,"w")
        f1.write(self.title)
        f1.write('1.0\n')
        for i in range(3) :
            for j in range(3) :
                f1.write(f"{self.lc[i][j]:23.16f}")
            f1.write('\n')
        for i in range(self.Ntype) :
            f1.write(self.atomtype[i]+"   ")
        f1.write('\n')
        for i in range(self.Ntype) :
            f1.write(str(self.Naint[i])+"  ")
        f1.write('\n')
        if self.f_seldyn :
           f1.write('Selective dynamics\n')
        f1.write('Direct\n')
        for i in range (self.Natom) :
            for j in range(3) :
                f1.write(f"{self.ap[i][j]:20.16f}")
            if self.f_seldyn :
                for j in range(3) :
                    if self.seldyn[i][j] :
                        f1.write(" T")
                    else :
                        f1.write(" F")
            f1.write('\n')
        f1.close()

    def match(self, Pright):
        for i in range(3):
            for j in range(3):
                if abs(self.lc[i][j]-Pright.lc[i][j])>0.0000001 :
                    print("Lattice vectors don't match.")
                    return False
        if self.Natom!=Pright.Natom or self.Ntype!=Pright.Ntype :
            print("Atoms don't match: #atoms of #types")
            return False
        for i in range(self.Ntype) :
            if self.atomtype[i]!=Pright.atomtype[i] or self.Naint[i]!=Pright.Naint[i] :
                print("Atoms don't match: "+str(i+1)+": "+self.atomtype[i]+" "+str(self.Naint[i])+" vs "+Pright.atomtype[i]+" "+str(Pright.Naint[i]))
                return False
        print("Structures match each other.")
        return True

    def moveatoms(self, Pright, factor):
        # move atoms by a new structure
        # Pright: POSCAR object
        for i in range(self.Natom) :
            for j in range(3) :
                if Pright.ap[i][j]-self.ap[i][j]>=0.5 :
                    newap=Pright.ap[i][j]-1.0
                elif Pright.ap[i][j]-self.ap[i][j]< -0.5 :
                    newap=Pright.ap[i][j]+1.0
                else :
                    newap=Pright.ap[i][j]
                self.ap[i][j]=self.ap[i][j]*(1.0-factor)+newap*factor
        self.movetobox()

    def movebyvector(self, disp, factor):
        # move atoms by a set of vectors
        # disp[Na][3]
        if len(disp)!=self.Natom:
            print("Warning: the length of displacement vector != Natom")
        fdisp=0.0
        for ia in range(self.Natom):
            fdisp += (disp[ia][0]**2+disp[ia][1]**2+disp[ia][2]**2)**0.5
        # the sum of all vectors are set to 1 A
        for ia in range(self.Natom):
            disp[ia][0]=disp[ia][0]/fdisp
            disp[ia][1]=disp[ia][1]/fdisp
            disp[ia][2]=disp[ia][2]/fdisp
        
        reclc=self.reclc_out()
        disp=np.matmul(np.array(disp),np.array(reclc)).tolist()
        print(disp[0])
        for ia in range(self.Natom) :
            for ix in range(3) :
                self.ap[ia][ix]=self.ap[ia][ix]+disp[ia][ix]*factor
        self.movetobox()

    def vacuum(self,newz) :
        factor=(newz/self.lc[2][2])
        if abs(self.lc[2][0])>0.000001 or abs(self.lc[2][1])>0.000001 :
            print("a3 doesn't parllel to z")
            sys.exit()
        self.movetobox()
        for i in range(self.Natom) :
            if self.ap[i][2]<=0.5 :
                self.ap[i][2]=self.ap[i][2]/factor
            else :
                self.ap[i][2]=1.0-(1.0-self.ap[i][2])/factor
        self.movetobox()
        self.lc[2][2]=newz

    def displacement(self, Pright):
        # find the displacement vector of two POSCARs
        # designed for a structural difference (poscar1.displacement(poscar2))
        disp=[]
        # disp[Na][3]
        for i in range(self.Natom) :
            disp0=[0.0,0.0,0.0]
            for j in range(3) :
                d=self.ap[i][j]-Pright.ap[i][j]
                while (True):
                   if d > 0.5+1.e-6 :
                       d=d-1.0
                   elif d < -0.5-1.e-6 :
                       d=d+1.0
                   else :
                       break
                disp0[0]=disp0[0]+d*self.lc[0][j]
                disp0[1]=disp0[1]+d*self.lc[1][j]
                disp0[2]=disp0[2]+d*self.lc[2][j]
            disp.append(disp0)
        return disp 

    def total_distance(self, Pright):
        # calculate the total distance (to zero) of a poscar
        # designed for a structural difference (poscar1.total_distance(poscar2))
        tot=0.0
        for i in range(self.Natom) :
            dis_i=[0.0,0.0,0.0]
            for j in range(3) :
                newap=self.ap[i][j]-Pright.ap[i][j]
                while (True):
                   if newap > 0.5+1.e-6 :
                       newap=newap-1.0
                   elif newap < -0.5-1.e-6 :
                       newap=newap+1.0
                   else :
                       break
                dis_i[0]=dis_i[0]+newap*self.lc[0][j]
                dis_i[1]=dis_i[1]+newap*self.lc[1][j]
                dis_i[2]=dis_i[2]+newap*self.lc[2][j]
            tot=tot+dis_i[0]**2+dis_i[1]**2+dis_i[2]**2
        tot=tot**0.5
        return tot 
            
    def general_q(self, Pright):
        # calculate the generalized cordinate of a total distance (to zero) of a poscar
        # designed for a structural difference (poscar1.general_q(poscar2))
        tot=0.0
        count=0
        it=0
        for i in range(self.Natom) :
            dis_i=[0.0,0.0,0.0]
            count+=1
            if count>self.Naint[it]:
                count=0
                it+=1
            for j in range(3) :
                newap=self.ap[i][j]-Pright.ap[i][j]
                while (True):
                   if newap > 0.5+1.e-6 :
                       newap=newap-1.0
                   elif newap < -0.5-1.e-6 :
                       newap=newap+1.0
                   else :
                       break
                dis_i[0]=dis_i[0]+newap*self.lc[0][j]
                dis_i[1]=dis_i[1]+newap*self.lc[1][j]
                dis_i[2]=dis_i[2]+newap*self.lc[2][j]
            tot=tot+(dis_i[0]**2+dis_i[1]**2+dis_i[2]**2)*self.dmass[self.atomtype[it]]
        tot=tot**0.5
        return tot 

    def new_seldyn(self) :
        self.f_seldyn=True
        self.seldyn=[]
        for i in range(self.Natom) :
            self.seldyn.append([True,True,True])

    def del_seldyn(self) :
        self.f_seldyn=False
        del self.seldyn

    def ap_out(self) :
        for i in range(self.Natom) :
            print(str(self.ap[i][0])+"  "+str(self.ap[i][1])+"  "+str(self.ap[i][2]))
        return(self.ap)

    def lc_out(self) :
        return(self.lc)

    def reclc_out(self) :
        reclc=np.linalg.inv(np.array(self.lc)).tolist()
        return(reclc)
   
    def atomtype_out(self) :
        return(self.atomtype)

    def Naint_out(self) :
        return(self.Naint)
