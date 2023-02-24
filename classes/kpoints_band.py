import os
import sys

class KPOINTS_band :
    def __init__(self, filename="KPOINTS",empty=False) :
        
        self.kpdict={} # a dictionary which maps label to b-coordinate
        self.kph=[] # a path of high symmetry k-points. use kpdict[kph[ik]][ix]

        if empty:
            self.nk_line=0
            return

        f0=open(filename,"r")
        line=f0.readlines()
        f0.close()
        if line[2][0] not in {"L","l"} :
            print("line-mode required.")
            sys.exit()
        if line[3][0] in {"C","c","K","k"} :
            print("fractional coordinates of k-points required.")
            sys.exit()

        self.nk_line=int(line[1].split()[0]) # kpoints on a line

        greek=self.loadgreek()

        for il in range(4,len(line)):
            word=line[il].split()
            if len(word)==0  or word[0][0]=="#" or word[0][0]=="!" :
                continue
            word=line[il].replace("!"," ").replace("#"," ").split()
            label=word[3]
            if label in greek:
                label=greek[label]
            if label not in self.kpdict :
                self.kpdict[label]=[float(word[0]),float(word[1]),float(word[2])]
            if len(self.kph)==0 or label != self.kph[-1] :
                self.kph.append(label)

    def fileread_kpathin(self,filename="kpath.in",nk_line=40):

        self.nk_line=nk_line

        f0=open(filename,"r")
        line=f0.readlines()
        f0.close()

        greek=self.loadgreek()

        for l in line:
            word=l.split()
            if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
                continue
            if len(word)==4 :
                label=word[3]
                if label in greek :
                    label=greek[label]
                self.kpdict[label]=[float(word[0]),float(word[1]),float(word[2])]
            elif len(word)==1 :
                label=word[0]
                if label in greek :
                    label=greek[label]
                self.kph.append(label)
            else :
                print("Error: check kpath.in")
                sys.exit()

#########################################################################

    def filewrite_vasp(self,filename="KPOINTS.new") :
        f1=open(filename,"w")
        f1.write("k points along high symmetry lines\n")
        f1.write("40\n")
        f1.write("line mode\n")
        f1.write("fractional\n")
        for ih in range(len(self.kph)-1) :
            for ix in range(3) :
                f1.write("{:16.12f}".format(self.kpdict[self.kph[ih]][ix])+" ")
            f1.write(self.kph[ih]+"\n")
            for ix in range(3) :
                f1.write("{:16.12f}".format(self.kpdict[self.kph[ih+1]][ix])+" ")
            f1.write(self.kph[ih+1]+"\n\n")
        f1.close()

    def filewrite_qe(self,filename="kpath.out"):
        kplot=[]
        kplot.append(self.kpdict[self.kph[0]])
        for ih in range(len(self.kph)-1):
            for ik in range(self.nk_line) :
                k=[0.0,0.0,0.0]
                for ix in range(3) :
                    k[ix]=self.kpdict[self.kph[ih]][ix]*float(self.nk_line-ik-1)/float(self.nk_line)+self.kpdict[self.kph[ih+1]][ix]*float(ik+1)/float(self.nk_line)
                kplot.append(k)

        f1=open(filename,"w")
        f1.write("K_POINTS crystal_b\n")
        f1.write(str(len(kplot))+"\n")
        for ik in range(len(kplot)) :
            for ix in range(3) :
                f1.write("{:16.12f}".format(kplot[ik][ix])+" ")
            f1.write("1.0\n")
        f1.close()

    def filewrite_kpathin(self, filename="kpath.in.new") :
        f1=open(filename,"w")
        for h in self.kpdict :
            for ix in range(3) :
                f1.write(str(self.kpdict[h][ix])+" ")
            f1.write(h+"\n")
        f1.write("\n")
        for h in self.kph :
            f1.write(h+"\n")
        f1.close()

    def filewrite_wannier90(self, filename="wannier90_kpath.dat") :
        f1=open(filename,"w")
        f1.write("begin kpoint_path\n")
        for ih in range(len(self.kph)-1) :
            f1.write(self.kph[ih]+" ")
            for ix in range(3) :
                f1.write("{:12.8f}".format(self.kpdict[self.kph[ih]][ix])+" ")
            f1.write(" "+self.kph[ih+1]+" ")
            for ix in range(3) :
                f1.write("{:12.8f}".format(self.kpdict[self.kph[ih+1]][ix])+" ")
            f1.write("\n")
        f1.write("end kpoint_path\n")
        f1.close()
            
#########################################################################

    def kph_out(self) :
        return(self.kph)

    def kphout_out(self,reclc):
        kphout=[]
        kphout.append(0.0)
        for ik in range(len(self.kph)) :
            if ik>0 :
                kpcold=[]
                kpcold.append(kpc[0])
                kpcold.append(kpc[1])
                kpcold.append(kpc[2])
                del kpc
            kpc=[]
            for i in range(3) :
                kpc0=0.0
                for j in range(3) :
                    kpc0=kpc0+self.kpdict[self.kph[ik]][j]*reclc[i][j]
                kpc.append(kpc0)
            if ik>0 :
                dkpc=((kpc[0]-kpcold[0])**2+(kpc[1]-kpcold[1])**2+(kpc[2]-kpcold[2])**2)**0.5
                kphout.append(kphout[ik-1]+dkpc)
                del kpcold
        return(kphout)

    def loadgreek(self):
        this_dir, this_filename = os.path.split(__file__)
        DATA_PATH = os.path.join(this_dir, "..", "data", "greek.dat")
        f0=open(DATA_PATH,"r")
        line=f0.readlines()
        f0.close()
        greek={}
        for l in line :
            word=l.split()
            if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
                continue
            if len(word) != 2:
                print("Error: check the structure of dftscr/data/greek.dat")
                sys.exit()
            greek[word[0]]=word[1]
        return greek
