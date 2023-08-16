import numpy as np
import sys
import math
import os
import xml.etree.ElementTree as ET

class POSCAR:
    # title
    # lc[3][3]
    # Natom
    # Ntype
    # f_seldyn
    # atomtype[Ntype]
    # Naint[Ntype]
    # ap[Natom][3]
    # seldyn[Natom][3] (if f_seldyn)
    # dmass{}
    
    def __init__(self, filename="POSCAR", empty=False) :
        if empty:
            self.title="SYSTEM\n"
            self.lc=[]
            self.Natom=0
            self.Ntype=0
            self.f_seldyn=False
            self.atomtype=[]
            self.Naint=[]
            self.ap=[]
            self.dmass={}
            return

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

    def fileread_qe(self, filename):
        f1=open(filename,"r")
        line=f1.readlines()
        f1.close()

        Fat=False

        for l in range(len(line)):
            word=line[l].split()
            if len(word)>=2 and word[0]=="CELL_PARAMETERS" :
                word[1]=word[1].replace("(","").replace("{","")
                if word[1][0]=="a" or word[1][0]=="A" :
                    factor=1.0
                elif word[1][0]=="b" or word[1][0]=="B" :
                    factor=0.529177249
                else :
                    print(word[1])
                    print("Only Angstrom and Bohr are supported.")
                    sys.exit()
                for il in range(3) :
                    wordlc=line[l+il+1].split()
                    self.lc.append([float(wordlc[0]),float(wordlc[1]),float(wordlc[2])])
                continue
            elif len(word)>=2 and word[0]=="ATOMIC_POSITIONS" :
                Fat=True
                continue
        
            if Fat==True :
                try :
                    self.ap.append([float(word[1]),float(word[2]),float(word[3])])
                except :
                    Fat=False
                    continue

                self.Natom+=1
                if self.Ntype==0 or self.atomtype[-1]!=word[0] :
                    self.Ntype+=1
                    self.atomtype.append(word[0])
                    self.Naint.append(1)
                else :
                    self.Naint[-1]+=1
        
        for i in range(3) :
            for j in range(3) :
                self.lc[i][j]=self.lc[i][j]*factor
        
        self.movetobox()

    def fileread_xml(self, filename=""):
        if filename=="" :
            # find a .xml file
            files = os.listdir()
            for f in files:
                if f.endswith('.xml'):
                    filename=f
                    break
        if filename=="" :
            print("Error: .xml file is not found")
            sys.exit()


        tree=ET.parse(filename)
        cell=tree.getroot().find("output").find("atomic_structure").find("cell")

        for ix in range(3):
            self.lc.append([float(lc1) for lc1 in cell.find("a"+str(ix+1)).text.split()])
        atoms=tree.getroot().find("output").find("atomic_structure").find("atomic_positions").findall("atom")
        
        for atom in atoms :
            atom_name=atom.get("name")
            if len(self.atomtype)==0 or self.atomtype[-1]!=atom_name:
                self.atomtype.append(atom_name)
                self.Naint.append(1)
            else :
                self.Naint[-1]+=1
            self.ap.append([float(ap) for ap in atom.text.split()])

        # convert Cartesian to direct
        self.ap=(np.array(self.ap)@np.linalg.inv(np.array(self.lc))).tolist()
        for ix in range(3):
            for iy in range(3):
                self.lc[ix][iy]=self.lc[ix][iy]*0.529177249

        self.Natom=len(self.ap)
        self.Ntype=len(self.atomtype)

    def fileread_prt(self, filename):
        f1=open(filename,"r")
        line=f1.readlines()
        f1.close()

        Flc=False
        Fat=False

        for l in range(len(line)):
            word=line[l].split()
            if len(word)>=2 and word[0]=="begin" and word[1]=="latticevecs" :
                Flc=True
                continue
            elif len(word)>=2 and word[0]=="end" and word[1]=="latticevecs" :
                Flc=False
                continue
            elif len(word)>=2 and word[0]=="begin" and word[1]=="coordinates" :
                Fat=True
                continue
            elif len(word)>=2 and word[0]=="end" and word[1]=="coordinates" :
                Fat=False
                continue
        
            if Flc==True :
                if len(word)>=4 and word[0]=="coord" :
                    self.lc.append([float(word[1]),float(word[2]),float(word[3])])
                elif len(word)>=2 and word[0]=="volume" :
                    volume=float(word[1])
            elif Fat==True :
                if len(word)>=2 and word[0]=="newtype" :
                    self.Ntype+=1
                    self.atomtype.append(word[1])
                    self.Naint.append(0)
                elif len(word)>=4 and word[0]=="coord" :
                    self.ap.append([float(word[1]),float(word[2]),float(word[3])])
                    self.Naint[-1]+=1
                    self.Natom+=1
        
        detlc=self.volume()
        factor=(volume/6.748334503468005/detlc)**(1.0/3.0)
        
        for i in range(3) :
            for j in range(3) :
                self.lc[i][j]=self.lc[i][j]*factor
        
        self.movetobox()

    def fileread_parsec(self, filename="parsec.in"):
        f1=open(filename,"r")
        line=f1.readlines()
        f1.close()

        Fat=False

        for il in range(len(line)):
            word=line[il].replace(":"," ").replace("="," ").split()
            if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
                continue
            if len(word)>=2 and word[0]=="begin" and word[1]=="atom_coord" :
                Fat=True
                continue
            if len(word)>=2 and word[0]=="end" and word[1]=="atom_coord" :
                Fat=False
                continue
            if word[0]=="Boundary_Sphere_Radius": 
                alat=float(word[1])*0.529177249*2.0 # in Bohr here
                self.lc=[[alat,0.0,0.0],[0.0,alat,0.0],[0.0,0.0,alat]]
            if word[0]=="atom_types_num" :
                self.Ntype=int(word[1])
            if word[0]=="atom_type" :
                self.atomtype.append(word[1])
                self.Naint.append(0)
            if word[0]=="coordinate_unit" and word[1]=="cartesian_ang" :
                factor=1.0
            if word[0]=="coordinate_unit" and word[1]=="cartesian_bohr" :
                factor=0.529177249
            if Fat==True :
                ap=[]
                for ix in range(3):
                    x=float(word[ix])*factor/alat+0.5
                    ap.append(x)
                self.ap.append(ap)
                self.Naint[-1]+=1
                self.Natom+=1

    def fileread_xyz(self,filename="") :
        if filename=="" :
            # find a .xyz file
            files = os.listdir()
            for f in files:
                if f.endswith('.xyz'):
                    filename=f
                    break
        if filename=="" :
            print("Error: .xyz file is not found")
            sys.exit()
        f0=open(filename,"r")
        line=f0.readlines()
        f0.close()

        radius=0.0
        self.Natom=int(line[0].split()[0])
        if len(line)<self.Natom+2 :
            print("Error: Length of file "+filename+" is less than Natom+2 "+str(self.Natom))
            sys.exit()
        for il in range(2,self.Natom+2) :
            word=line[il].split()
            if self.Ntype==0 or self.atomtype[-1]!=word[0] :
                self.Ntype+=1
                self.atomtype.append(word[0])
                self.Naint.append(1)
            else :
                self.Naint[-1]+=1
            ap0=[float(word[1]),float(word[2]),float(word[3])]
            self.ap.append(ap0)
            radius=max(radius,abs(ap0[0]),abs(ap0[1]),abs(ap0[2]))
            
        lc0=2.0*radius+10.0
        self.lc=[[lc0,0,0],[0,lc0,0],[0,0,lc0]]
        for ia in range(self.Natom) :
            for ix in range(3) :
                self.ap[ia][ix]=self.ap[ia][ix]/lc0+0.5
            
########################################################################

    def filewrite(self, filename="POSCAR.new"):
        f1=open(filename,"w")
        f1.write(self.title)
        f1.write('1.0\n')
        for i in range(3) :
            for j in range(3) :
                f1.write(f"{self.lc[i][j]:21.12f}")
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
                f1.write(f"{self.ap[i][j]:18.12f}")
            if self.f_seldyn :
                for j in range(3) :
                    if self.seldyn[i][j] :
                        f1.write(" T")
                    else :
                        f1.write(" F")
            f1.write('\n')
        f1.close()

    def filewrite_qe(self, filename="qe_st.dat"):
        self.load_dmass()

        kgrid=[]
        for i in range(3) :
            kgrid.append(math.floor(30.0/(self.lc[i][0]**2+self.lc[i][1]**2+self.lc[i][2]**2)**0.5)+1)
   
        f2=open(filename,"w")
        f2.write("CELL_PARAMETERS angstrom\n")
        for i in range(3) :
            f2.write(f"  {self.lc[i][0]:.12f}  {self.lc[i][1]:.12f}  {self.lc[i][2]:.12f}\n")
        f2.write("ATOMIC_SPECIES\n")
        for i in range(self.Ntype) :
            f2.write("  "+self.atomtype[i]+"  "+str(self.dmass[self.atomtype[i]])+"  "+self.atomtype[i]+".upf\n")
        f2.write("ATOMIC_POSITIONS crystal\n")
        ij=0
        ik=0
        for i in range(self.Natom) :
            f2.write(f"  {self.atomtype[ij]:2s}  {self.ap[i][0]:.12f}  {self.ap[i][1]:.12f}  {self.ap[i][2]:.12f}\n")
            ik=ik+1
            if ik==self.Naint[ij] :
                ij=ij+1
                ik=0
        f2.write("K_POINTS automatic\n")
        f2.write(f"  {kgrid[0]:d}  {kgrid[1]:d}  {kgrid[2]:d} 0 0 0\n")
        f2.close()

    def filewrite_prt(self, filename="prt_st.dat"):
        volume=self.volume()*6.748334503468005
        # convert Angstrom^3 to bohr^3
        
        f2=open(filename,"w")
        f2.write("begin latticevecs\n")
        for i in range(3) :
            f2.write("coord  "+f"{self.lc[i][0]:20.12f}"+"  "+f"{self.lc[i][1]:20.12f}"+"  "+f"{self.lc[i][2]:20.12f}"+"\n")
        f2.write("volume "+f"{volume:24.16f}"+"\nend latticevecs\n\nbegin coordinates\n")
        k=0
        for i in range(self.Ntype) :
            f2.write("newtype "+self.atomtype[i]+"\n")
            for j in range(self.Naint[i]) :
                f2.write("coord  "+f"{self.ap[k][0]:20.12f}"+"  "+f"{self.ap[k][1]:20.12f}"+"  "+f"{self.ap[k][2]:20.12f}"+"\n")
                k=k+1
        f2.write("end coordinates\n\n")
        f2.close()

    def filewrite_parsec(self, filename="parsec_st.dat"):
        f2=open(filename,"w")
        f2.write("begin cell_shape\n")
        for i in range(3) :
            f2.write("  "+f"{self.lc[i][0]:20.12f}"+"  "+f"{self.lc[i][1]:20.12f}"+"  "+f"{self.lc[i][2]:20.12f}"+"\n")
        f2.write("end cell_shape\n\n")
        k=0
        for i in range(self.Ntype) :
            f2.write("#------------- new atom type -------------\n")
            f2.write("Atom_Type: "+self.atomtype[i]+"\n")
            #f2.write("Pseudopotential_Format: \n")
            #f2.write("Core_Cutoff_Radius: \n")
            f2.write("Local_Component: \n")
            #f2.write("Potential_Num: \n\n")
            #f2.write("begin Electron_Per_Orbital\n")
            #f2.write("# S P D F\n \n")
            #f2.write("end Electron_Per_Orbital\n\n")

            f2.write("begin Atom_Coord\n")
            for j in range(self.Naint[i]) :
                f2.write("  "+f"{self.ap[k][0]:20.12f}"+"  "+f"{self.ap[k][1]:20.12f}"+"  "+f"{self.ap[k][2]:20.12f}"+"\n")
                k=k+1
            f2.write("end Atom_Coord\n\n")
        f2.close()

    def filewrite_wannier90(self, filename="wannier90_st.dat"):
        f2=open(filename,"w")
        f2.write("Begin Unit_Cell_Cart\n")
        for i in range(3) :
            f2.write(f"  {self.lc[i][0]:.12f}  {self.lc[i][1]:.12f}  {self.lc[i][2]:.12f}\n")
        f2.write("End Unit_Cell_Cart\n\n")

        f2.write("Begin Projections\n")
        f2.write("  random\n")
        for i in range(self.Ntype) :
            f2.write("  "+self.atomtype[i]+"  :\n")
        f2.write("End Projections\n\n")

        f2.write("Begin Atoms_Frac\n")
        ij=0
        ik=0
        for i in range(self.Natom) :
            f2.write(f"  {self.atomtype[ij]:2s}  {self.ap[i][0]:.12f}  {self.ap[i][1]:.12f}  {self.ap[i][2]:.12f}\n")
            ik=ik+1
            if ik==self.Naint[ij] :
                ij=ij+1
                ik=0
        f2.write("End Atoms_Frac\n\n")
        f2.close()

    def filewrite_xyz(self, filename="structure.xyz") :
        f2=open(filename,"w")
        f2.write(str(self.Natom)+"\n\n")
        it1=0
        it2=0
        for ia in range(self.Natom) :
            apc=[0.0]*3
            for ix in range(3) :
                apc[0]+=(self.ap[ia][ix]-0.5)*self.lc[ix][0]
                apc[1]+=(self.ap[ia][ix]-0.5)*self.lc[ix][1]
                apc[2]+=(self.ap[ia][ix]-0.5)*self.lc[ix][2]
            f2.write(f"  {self.atomtype[it1]:2s}{apc[0]:16.10f}{apc[1]:16.10f}{apc[2]:16.10f}\n")
            it2=it2+1
            if it2==self.Naint[it1] :
                it1=it1+1
                it2=0
            
        f2.close()

########################################################################

    def movetobox(self) :
        for i in range(self.Natom) :
            for j in range(3) :
                self.ap[i][j]=self.ap[i][j]-self.ap[i][j]//1.0
                if self.ap[i][j]<0.0000000001 or self.ap[i][j]>0.9999999999 :
                    self.ap[i][j]=0.0

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
        # The vector is defined in a cartesian coordinate system using angstrom.
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
        self.load_dmass()
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


    def reclc_out(self) :
        reclc=np.linalg.inv(np.array(self.lc)).tolist()
        for i in range(3) :
            for j in range(3) :
                reclc[i][j]=reclc[i][j]*2.0*3.141592653589793238463 
        return(reclc)
   
    def volume(self) :
        return  self.lc[0][0]*(self.lc[1][1]*self.lc[2][2]-self.lc[1][2]*self.lc[2][1])+self.lc[0][1]*(self.lc[1][2]*self.lc[2][0]-self.lc[1][0]*self.lc[2][2])+self.lc[0][2]*(self.lc[1][0]*self.lc[2][1]-self.lc[1][1]*self.lc[2][0])

    def load_dmass(self) :
        if self.dmass :
            return
        this_dir, this_filename = os.path.split(__file__)
        DATA_PATH = os.path.join(this_dir, "..", "data", "atomicmass.dat")
        f0=open(DATA_PATH,"r")
        line=f0.readlines()
        f0.close()

        for l in line:
            word=l.split()
            if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
                continue
            self.dmass[word[2]]=float(word[3])
