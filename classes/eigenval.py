import numpy as np
import xml.etree.ElementTree as ET

class EIGENVAL : 
    # Nk, Nb, Ns
    # eig[k][b][s]
    # occ[k][b][s]
    # Nvb[s]
    # is_semic
    # kp[k][3]
    # wt[k][3]

    def __init__(self, filename="EIGENVAL",is_hse=False,empty=False) :
        self.wt=[]
        self.kp=[]
        self.eig=[]
        self.occ=[]

        if empty:
            self.Nk=0
            self.Nb=0
            self.Ns=1
            return

        f0=open(filename,"r")
        line=f0.readlines()
        f0.close()
        self.Ns=int(line[0].split()[3])
        self.Nk=int(line[5].split()[1])
        self.Nb=int(line[5].split()[2])
        for k in range(self.Nk) :
            word=line[7+k*(self.Nb+2)].split()
            self.kp.append([float(word[0]),float(word[1]),float(word[2])])
            self.wt.append(float(word[3]))
            eig0=[]
            occ0=[]
            for b in range(self.Nb) :
                word=line[8+k*(self.Nb+2)+b].split()
                if self.Ns==1 :
                    eig0.append([float(word[1])])
                    occ0.append([float(word[2])])
                elif self.Ns==2 :
                    eig0.append([float(word[1]),float(word[2])])
                    occ0.append([float(word[3]),float(word[4])])
            self.eig.append(eig0)
            self.occ.append(occ0)

            if (is_hse and abs(self.wt[len(self.wt)-1]) > 0.0000001) :
                self.kp.pop(len(self.kp)-1)
                self.wt.pop(len(self.wt)-1)
                self.eig.pop(len(self.eig)-1)
                self.occ.pop(len(self.occ)-1)

        if (is_hse) :
            self.Nk=len(self.kp)

        self.semic()
        if self.is_semic==True :
            self.gap()
                    
        del line

    def fileread_qexml(self, filename="*.xml") :
        # only support Ns=1 here
        # Nk, Nb, Ns
        # eig[k][b][s]
        # occ[k][b][s]
        # Nvb[s]
        # is_semic
        # kp[k][3]
        # wt[k]
        tree=ET.parse(filename)
        kpoints=tree.getroot().find("output").find("band_structure").findall("ks_energies")

        self.Nk=len(kpoints)
        self.Nb=int(kpoints[0].find('eigenvalues').get("size"))
        self.Ns=1
        
        for kp in kpoints :
            # use the unit eV and A (^-1)
            self.kp.append([float(kp1)/0.529177210903 for kp1 in kp.find('k_point').text.split()])
            self.wt.append(float(kp.find('k_point').get("weight")))
            # qe xml energies are in Hatree
            self.eig.append([[float(eig1)*27.211386245988] for eig1 in kp.find('eigenvalues').text.split()])
            self.occ.append([[float(occ1)] for occ1 in kp.find('occupations').text.split()])

        self.semic()
        if self.is_semic==True :
            self.gap()

#########################################################################
    
    def eigshift(self,ezero) :
        for k in range(self.Nk) :
            for b in range(self.Nb) :
                for s in range(self.Ns) :
                    self.eig[k][b][s]=self.eig[k][b][s]-ezero

    def bandkpout(self,reclc=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]) :
        # transfer the fractional to cartesian in k space
        # if kp is already in cartisian, simply set reclc=I
        kpout=[]
        kpout.append(0.0)
        for k in range(self.Nk) :
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
                    kpc0=kpc0+self.kp[k][j]*reclc[i][j]
                kpc.append(kpc0)
            if k>0 :
                dkpc=((kpc[0]-kpcold[0])**2+(kpc[1]-kpcold[1])**2+(kpc[2]-kpcold[2])**2)**0.5
                kpout.append(kpout[k-1]+dkpc)
                del kpcold
        return(kpout)

    def eigtrans(self) :
        eigout=[]
        for s in range(self.Ns) :
            eigout1=[]
            for b in range(self.Nb) :
                eigout0=[]
                for k in range(self.Nk) :
                    eigout0.append(self.eig[k][b][s])
                eigout1.append(eigout0)
            eigout.append(eigout1)
        return(eigout)

    def writegap(self, kp) :
        f0=open("gap.txt","w")
        if self.is_semic==True :
            if self.vbm_k!=self.cbm_k :
                if self.vbm_k%kp.nk_line==0 :
                    vbm_kl=kp.kph[self.vbm_k//kp.nk_line]
                else :
                    vbm_kl="("+kp.kph[self.vbm_k//kp.nk_line]+","+kp.kph[self.vbm_k//kp.nk_line+1]+")"
                if self.cbm_k%kp.nk_line==0 :
                    cbm_kl=kp.kph[self.cbm_k//kp.nk_line]
                else :
                    cbm_kl="("+kp.kph[self.cbm_k//kp.nk_line]+","+kp.kph[self.cbm_k//kp.nk_line+1]+")"
                eindg_print="Indirect: Eg = "+str(round(self.eindg,4))+" eV, between "+vbm_kl+" -> "+cbm_kl
                f0.write(eindg_print+"\n")
                print(eindg_print)
            if self.edg_k%kp.nk_line==0 :
                edg_kl=kp.kph[self.edg_k//kp.nk_line]
            else :
                edg_kl="("+kp.kph[self.edg_k//kp.nk_line]+","+kp.kph[self.edg_k//kp.nk_line+1]+")"
            edg_print="Direct:   Eg = "+str(round(self.edg,4))+" eV, at "+edg_kl
            f0.write(edg_print+"\n")
            print(edg_print)
        else :
            f0.write("Eg = 0, no band gap\n")
        f0.close()

    def semic(self) :
        self.is_semic=True
        for k in range(self.Nk) :
            for b in range(self.Nb) :
                for s in range(self.Ns) :
                    if self.occ[k][b][s]<0.999 and self.occ[k][b][s]>0.001 :
                        self.is_semic=False
    def gap(self) :
        self.Nvb=[]
        for s in range(self.Ns) :
            self.Nvb.append(0)
            for b in range(self.Nb) :
                if self.occ[0][b][s]<0.5 :
                    self.Nvb[s]=b
                    break
        self.vbm=-1e10
        self.vbm_k=0
        self.vbm_s=0
        self.cbm=1e10
        self.cbm_k=0
        self.cbm_s=0
        self.edg=1e10
        self.edg_k=0
        self.edg_s=0
        for s in range(self.Ns) :
            for k in range(self.Nk) :
                if self.eig[k][self.Nvb[s]-1][s]>self.vbm :
                    self.vbm=self.eig[k][self.Nvb[s]-1][s]
                    self.vbm_k=k
                    self.vbm_s=s
                if self.eig[k][self.Nvb[s]][s]<self.cbm :
                    self.cbm=self.eig[k][self.Nvb[s]][s]
                    self.cbm_k=k
                    self.cbm_s=s
                if self.eig[k][self.Nvb[s]][s]-self.eig[k][self.Nvb[s]-1][s]<self.edg :
                    self.edg=self.eig[k][self.Nvb[s]][s]-self.eig[k][self.Nvb[s]-1][s]
                    self.edg_k=k
                    self.edg_s=s
        self.eindg=self.cbm-self.vbm