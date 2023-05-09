import os
import sys

class KPOINTS_band :
    def __init__(self, filename="KPOINTS",empty=False) :
        
        self.kpdict={} # a dictionary which maps label to b-coordinate
        self.kph=[] # a path of high symmetry k-points. kph[Np][Nh[ip]]
        # use kpdict[kph[ip][ik]][ix] for the ik-th points in the ip-th path. ix indicates kx, ky, or kz

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

        ip=0
        ffirst=True
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
            if ffirst==True:
                if len(self.kph)==0 or label != self.kph[-1][-1]:
                    self.kph.append([label])
            else:
                self.kph[-1].append(label)
            ffirst=not ffirst

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
            #if len(word)==4 :
            try :
                float(word[0])
                float(word[1])
                float(word[2])
            except :
                self.kph.append(word)
                for ik in range(len(self.kph[-1])) :
                    if self.kph[-1][ik] in greek :
                        self.kph[-1][ik]=greek[self.kph[-1][ik]]
            else :
                label=word[3]
                if label in greek :
                    label=greek[label]
                self.kpdict[label]=[float(word[0]),float(word[1]),float(word[2])]

#########################################################################

    def filewrite_vasp(self,filename="KPOINTS.new") :
        f1=open(filename,"w")
        f1.write("k points along high symmetry lines\n")
        f1.write("40\n")
        f1.write("line mode\n")
        f1.write("fractional\n")
        for p in self.kph :
            for ih in range(len(p)-1) :
                for ix in range(3) :
                    f1.write("{:16.12f}".format(self.kpdict[p[ih]][ix])+" ")
                f1.write(p[ih]+"\n")
                for ix in range(3) :
                    f1.write("{:16.12f}".format(self.kpdict[p[ih+1]][ix])+" ")
                f1.write(p[ih+1]+"\n\n")
        f1.close()

    def filewrite_qe(self,filename="kpath.out"):
        kplot=[]
        for ip in range(len(self.kph)):
            kplot.append(self.kpdict[self.kph[ip][0]])
            for ih in range(len(self.kph[ip])-1):
                for ik in range(self.nk_line) :
                    k=[0.0,0.0,0.0]
                    for ix in range(3) :
                        k[ix]=self.kpdict[self.kph[ip][ih]][ix]*float(self.nk_line-ik-1)/float(self.nk_line)+self.kpdict[self.kph[ip][ih+1]][ix]*float(ik+1)/float(self.nk_line)
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
        for p in self.kph :
            for h in p :
                f1.write(h+" ")
            f1.write("\n")
        f1.close()

    def filewrite_wannier90(self, filename="wannier90_kpath.dat") :
        f1=open(filename,"w")
        f1.write("bands_plot = true\n\n")
        f1.write("begin kpoint_path\n")
        for p in self.kph :
            for ih in range(len(p)-1) :
                f1.write(p[ih]+" ")
                for ix in range(3) :
                    f1.write("{:12.8f}".format(self.kpdict[p[ih]][ix])+" ")
                f1.write(" "+p[ih+1]+" ")
                for ix in range(3) :
                    f1.write("{:12.8f}".format(self.kpdict[p[ih+1]][ix])+" ")
                f1.write("\n")
        f1.write("end kpoint_path\n\n")
        f1.close()
            
#########################################################################

    def kphlabel_out(self):
        kphlabel=self.kph[0]
        for ip in range(1,len(self.kph)) :
            kphlabel[-1]=kphlabel[-1]+"|"+self.kph[ip][0]
            kphlabel=kphlabel+self.kph[ip][1:]
        return kphlabel

    def kphx_out(self,reclc):
        kphx=[0.0]
        for p in self.kph :
            for ik in range(len(p)) :
                if ik>0 :
                    kpcold=[kpc[0],kpc[1],kpc[2]]
                    del kpc
                kpc=[]
                for i in range(3) :
                    kpc0=0.0
                    for j in range(3) :
                        kpc0=kpc0+self.kpdict[p[ik]][j]*reclc[i][j]
                    kpc.append(kpc0)
                if ik>0 :
                    dkpc=((kpc[0]-kpcold[0])**2+(kpc[1]-kpcold[1])**2+(kpc[2]-kpcold[2])**2)**0.5
                    kphx.append(kphx[-1]+dkpc)
                    del kpcold
        return(kphx)

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

    def findlabel(self,k,dim=0) :
        # k[3] is a k point in direct coordinate,
        # this function find the label of k
        # dim is the maximum dimension of searching. (point, line, plane, ...)
        if dim>=0 :
            for h in self.kpdict :
                kh=self.kpdict[h]
                if abs(kh[0]-k[0])<=1e-5 and abs(kh[1]-k[1])<=1e-5 and abs(kh[2]-k[2])<=1e-5 :
                    return h
        if dim>=1 :
            for h1 in self.kpdict :
                kh1=self.kpdict[h1]
                for h2 in self.kpdict :
                    if h1==h2 :
                        continue
                    kh2=self.kpdict[h2]
                    x1=kh1[0]-k[0]
                    x2=kh2[0]-k[0]
                    y1=kh1[1]-k[1]
                    y2=kh2[1]-k[1]
                    z1=kh1[2]-k[2]
                    z2=kh2[2]-k[2]
                    if abs(x1*y2-x2*y1)<=1e-5 and abs(y1*z2-y2*z1)<=1e-5 :
                        return "("+h1+","+h2+")"
        return "elsewhere"

