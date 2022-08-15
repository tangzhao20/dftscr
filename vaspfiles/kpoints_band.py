import numpy as np

class KPOINTS_band :
    def __init__(self, filename="KPOINTS") :
        f0=open(filename,"r")
        line=f0.readlines()
        f0.close()
        if line[2][0]!="L" and line[2][0]!="l" :
            print("line-mode required.")
            sys.exit()
        if line[3][0]!="R" and line[3][0]!="r" :
            print("rec required.")
            sys.exit()

        self.nk_line=int(line[1].split()[0])
        Flaglinebegin=1
        self.kph=[]
        self.kphl=[]
        Flagfirstpoint=1
        Flagoutput=0
        for i in range(4,len(line)):
            word=line[i].replace("!"," ").split()
            if len(word)==0 :
                continue
            if Flagfirstpoint==1 :
                self.kph.append([float(word[0]),float(word[1]),float(word[2])])
                Flagoutput=1
                Flagfirstpoint=0
            if Flaglinebegin==0 : 
                self.kph.append([float(word[0]),float(word[1]),float(word[2])])
                Flagoutput=1
            if Flagoutput==1 :
                if word[3]=="G" or word[3]=="Gamma" :
                    self.kphl.append("Γ")
                else :
                    self.kphl.append(word[3])
                Flagoutput=0
            Flaglinebegin=1-Flaglinebegin


    def kphl_out(self) :
        return(self.kphl)

    def kphout_out(self,reclc):
        kphout=[]
        kphout.append(0.0)
        for k in range(len(self.kph)) :
            if k>0 :
                kpcold=[]
                kpcold.append(kpc[0])
                kpcold.append(kpc[1])
                kpcold.append(kpc[2])
                del kpc
            kpc=[]
            for i in range(3) :
                kpc0=0.0
                for j in range(3) :
                    kpc0=kpc0+self.kph[k][j]*reclc[i][j]
                kpc.append(kpc0)
            if k>0 :
                dkpc=((kpc[0]-kpcold[0])**2+(kpc[1]-kpcold[1])**2+(kpc[2]-kpcold[2])**2)**0.5
                kphout.append(kphout[k-1]+dkpc)
                del kpcold
        return(kphout)

