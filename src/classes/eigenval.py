import os
import numpy as np
import xml.etree.ElementTree as ET
from commons import load_constant

class EIGENVAL : 
    # Nk, Nb, Ns
    # eig[k][b][s]
    # occ[k][b][s]
    # Nvb[s]
    # is_semic
    # kp[k][3]
    # wt[k]

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

    def __str__(self) :
        str_out = "EIGENVAL:\n"
        str_out += " Nk = " + str(self.Nk) + "\n"
        str_out += " Nb = " + str(self.Nb) + "\n"
        str_out += " Ns = " + str(self.Ns) + "\n"
        if self.is_semic :
            str_out += " is_semic = True\n"
            str_out += f" vbm = {self.vbm:.2f} eV\n"
            str_out += f" cbm = {self.cbm:.2f} eV\n"
            str_out += f" Egdir = {self.edg:.2f} eV\n"
            if self.eindg < self.edg : 
                str_out += f" Egind = {self.eindg:.2f} eV\n"
        else :
            str_out += " is_semic = False\n"
        return str_out

#########################################################################

    def fileread_qexml(self, filename="") :
        Ha=load_constant("rydberg")*2.0

        if filename=="" :
            # find a .xml file
            files = os.listdir()
            for f in files:
                if f.endswith('.xml'):
                    filename=f
                    break
        tree=ET.parse(filename)
        kpoints=tree.getroot().find("output").find("band_structure").findall("ks_energies")
        cell=tree.getroot().find("input").find("atomic_structure").find("cell")
        lc=[]
        for ix in range(3):
            lc.append([float(lc1) for lc1 in cell.find("a"+str(ix+1)).text.split()])
            # convert lc (from Bohr) to a/2pi
        factor=(lc[0][0]**2+lc[0][1]**2+lc[0][2]**2)**0.5
        for ix1 in range(3) :
            for ix2 in range(3) :
                lc[ix1][ix2]=lc[ix1][ix2]/factor

        self.Nk=len(kpoints)
        if tree.getroot().find("output").find("band_structure").find("lsda").text in ["True", "true"] :
            self.Ns=2
        else :
            self.Ns=1
        self.Nb=int(kpoints[0].find('eigenvalues').get("size"))//self.Ns
        
        for kp in kpoints :
            # energies are in eV, and 2pi/a
            self.kp.append([float(kp1) for kp1 in kp.find('k_point').text.split()])
            self.wt.append(float(kp.find('k_point').get("weight")))
            eig1=kp.find('eigenvalues').text.split()
            occ1=kp.find('occupations').text.split()
            # qe xml energies are in Hatree (Ha)
            eig0=[]
            occ0=[]
            if self.Ns==1 :
                for ib in range(self.Nb) :
                    eig0.append([float(eig1[ib])*Ha])
                    occ0.append([float(occ1[ib])])
            else : # Ns==2
                for ib in range(self.Nb) :
                    eig0.append([float(eig1[ib])*Ha,float(eig1[ib+self.Nb])*Ha])
                    occ0.append([float(occ1[ib]),float(occ1[ib+self.Nb])])
            self.eig.append(eig0)
            self.occ.append(occ0)

        self.kc2kd(lc)

        self.semic()
        if self.is_semic==True :
            self.gap()

    def fileread_wan(self, Nb_pad=0) :
        # only support Ns=1 and semiconductor
        files = os.listdir()

        for f in files:
            if f.endswith('_band.kpt'):
                filename=f
                break
        f1=open(filename,"r")
        line=f1.readlines()
        f1.close()
        self.Nk=int(line[0].split()[0])
        for ik in range(self.Nk) :
            word=line[ik+1].split()
            kpoint=[]
            for ix in range(3) :
                kpoint.append(float(word[ix]))
            self.kp.append(kpoint)
            self.wt.append(float(word[3]))

        # find a *_band.dat file
        for f in files:
            if f.endswith('_band.dat'):
                filename=f
                break
        f2=open(filename,"r")
        line=f2.readlines()
        f2.close()

        for ik in range(self.Nk) :
            self.eig.append([])
        ik=0
        ib=-1
        fempty=True
        for l in line:
            word=l.split()
            if len(word)==2 :
                if fempty==True :
                    ik=0
                    ib=ib+1
                    fempty=False
                self.eig[ik].append([float(word[1])])
                ik=ik+1
            else :
                fempty=True
        self.Nb=len(self.eig[0])
        self.Ns=1

    def fileread_parsec(self) :
        bohr=load_constant("bohr")
        pi=load_constant("pi")
        rydberg=load_constant("rydberg")

        if os.path.isfile("bands.dat") :
            filename="bands.dat"
        elif os.path.isfile("eigen.dat"):
            self.fileread_parsec_eigen()
            return
        else :
            print("Error: Reading parsec eigenval needs bands.dat or eigen.dat")
            sys.exit()

        f1=open(filename,"r")
        line=f1.readlines()
        f1.close()

        word=line[1].split()
        self.Ns=int(word[0])
        if self.Ns==2 :
            print("Error: Ns==2 is not supported for parsec. See eigenval.py")
            sys.exit()
        self.Nb=int(word[1])
        Np=int(word[2])
        ef=float(word[3])
        self.Nk=0
        for ip in range(Np) :
            word=line[ip+3].split()
            self.Nk+=int(word[1])

        il=4+Np
        for ik in range(self.Nk): 
            word=line[il].split()
            self.kp.append([float(word[3]),float(word[4]),float(word[5])])
            il+=1
            eig0=[]
            for ib in range(self.Nb) :
                eig0.append([float(line[il].split()[0])])
                il+=1
            self.eig.append(eig0)

        # convert kp from 1/bohr to 1/A
        for ik in range(self.Nk) :
            for ix in range(3) :
                self.kp[ik][ix]=self.kp[ik][ix]/bohr*0.5/pi

        self.eigshift(ef)
        Nvb=[]
        for ik in range(self.Nk) :
            for ib in range(self.Nb) :
                if self.eig[ik][ib][0]>1e-6 :
                    Nvb.append(ib)
                    break
        if max(Nvb)==min(Nvb) :
            self.is_semic=True
            self.Nvb=[max(Nvb)]
        else :
            self.is_semic=False
        for ik in range(self.Nk) :
            for ib in range(self.Nb) :
                self.eig[ik][ib][0]*=rydberg

    def fileread_parsec_eigen(self) :
        # Read the eigenvalues from eigen.dat.
        # kp is missing in this file, so we assume the first k is gamma and read only this point.
        rydberg=load_constant("rydberg")

        filename="eigen.dat"

        f0=open(filename,"r")
        line=f0.readlines()
        f0.close()

        self.Nk=1
        self.Ns=0
        self.eig=[[]]
        self.occ=[[]]
        for l in line :
            word=l.split() 
            if len(word)==0 or word[0][0] in ["!","#"] :
                continue
            if not word[0].isdigit() :
                continue
            if self.Ns==0 :
                if len(word)==7 :
                    self.Ns=1
                elif len(word)==8 :
                    self.Ns=2
                else :
                    print("Error: number of columns in eigen.dat should be 7 or 8")
                    sys.exit()
            ik=int(word[5])-1
            if ik>0 :
                continue
            ib=int(word[0])-1
            eig0=float(word[1])*rydberg
            occ0=float(word[3])
            if self.Ns==2 and word[7]=="dn" :
                self.eig[0][ib].append(eig0)
                self.occ[0][ib].append(occ0)
            else : # Ns=1 or spin=up
                self.eig[0].append([eig0])
                self.occ[0].append([occ0])
        self.Nb=len(self.eig[0])
        self.kp=[[0.0,0.0,0.0]]
        self.wt=[[1.0]]
        self.semic()
        self.gap()

#########################################################################
    
    def eigshift(self,ezero) :
        for k in range(self.Nk) :
            for b in range(self.Nb) :
                for s in range(self.Ns) :
                    self.eig[k][b][s]=self.eig[k][b][s]-ezero

    def bandkpout(self,kp,reclc=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]) :
        # transfer the fractional to cartesian in k space
        # if kp is already in cartisian, simply set reclc=I
        kpc=np.array(self.kp)@(np.array(reclc).T)
        kplabelold=""
        kpout=[]
        for ik in range(self.Nk) :
            kplabel=kp.findlabel(self.kp[ik],dim=0)
            if kplabel!="elsewhere" and kplabel!=kplabelold and (kplabelold!="elsewhere" or max([abs(self.kp[ik][ix]-self.kp[ik-1][ix]) for ix in range(3)])>0.1) :
                kpout.append([0.0])
            else :
                dkpc=sum((kpc[ik]-kpc[ik-1])**2)**0.5
                kpout[-1].append(kpout[-1][-1]+dkpc)
            kplabelold=kplabel
        return(kpout)

    def kc2kd(self,lc) :
        #convert kc from Cartesian to direct
        self.kp=(np.array(self.kp)@(np.array(lc).T)).tolist()

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
                vbm_kl=kp.findlabel(self.kp[self.vbm_k],dim=1)
                cbm_kl=kp.findlabel(self.kp[self.cbm_k],dim=1)
                eindg_print="Indirect: Eg = "+str(round(self.eindg,4))+" eV, between "+vbm_kl+" -> "+cbm_kl
                f0.write(eindg_print+"\n")
                print(eindg_print)
            edg_kl=kp.findlabel(self.kp[self.edg_k],dim=1)
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
        # edg: direct band gap
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

