#!/usr/bin/env python3
# This tool prints some information from parsec wavefunction file parsec.dat
# Now only support the real wavefunction of cluster, which has 1 k-point (Gamma)

# python3 wfninfo.py

# outputs (wfninfo.out) are shown on the screen

from commons import load_constant
import sys
import os
import numpy as np
from v3math import v3tm3


def formatprint(mat, dim0, fstr) :
    # dim0 is the number of rows
    # example of fstr: "{:3}" or "{:.2f}"
    dim1=len(mat)//dim0
    strformat=dim1*fstr
    for ii in range(dim0) :
        f1.write(strformat.format(*mat[ii*dim1:(ii+1)*dim1]))
        f1.write("\n")

Ry=load_constant("rydberg")

f0=open("parsec.dat","rb")
f1=open("wfninfo.out","w")

datelabel=np.fromfile(f0, dtype='S30', count=1)[0][4:-5].decode()
f0.read(8)
f1.write(datelabel+"\n")

stype=np.fromfile(f0, dtype=np.int32, count=10)
f0.read(8)
f1.write("Nspin = "+str(stype[0])+"\n")
if stype[1]==0 :
    f1.write("Real wavefunctions\n")
elif stype[1]==1 :
    f1.write("Complex wavefunctions\n")
elif stype[1]==2 :
    f1.write("Spin-orbit (complex doubled wavefunctions)\n")
elif stype[1]==3 :
    f1.write("Noncollinear (spin3d present)\n")
elif stype[1]==4 :
    f1.write("Spin-orbit and noncollinear (spin3d present)\n")
if stype[2]==0 :
    f1.write("Cluster, no periodicity\n")
elif stype[2]==1 :
    f1.write("Bulk, periodicity along all cartesian directions\n")
elif stype[2]==2 :
    f1.write("Wire, periodicity along x direction only\n")
elif stype[2]==3 :
    f1.write("Slab, periodicity along x and y directions\n")

# pbc is the periodic boundary condition
if stype[2]==0 :
    # if cluster
    gridstep=np.fromfile(f0, dtype=np.double, count=1)[0]
    boundary=np.fromfile(f0, dtype=np.double, count=1)[0]
    f0.read(8)
    f1.write("Grid step = "+str(gridstep)+" Bohr\n")
    f1.write("Boundary size = "+str(boundary)+" Bohr\n")
else :
    # if pbc present
    print("Error: pbc present, wavinfo.py need work.")
    sys.exit()

Ndim=np.fromfile(f0, dtype=np.int32, count=1)[0]
f0.read(8)
f1.write("Ndim = "+str(Ndim)+"\n")

shift=np.fromfile(f0, dtype=np.double, count=3)
f0.read(8)

nwedge=np.fromfile(f0, dtype=np.int32, count=1)[0]
ntrans=np.fromfile(f0, dtype=np.int32, count=1)[0]
f0.read(8)
f1.write("Nwedge = "+str(nwedge)+" (Ndim in irreducible wedge)\n")
f1.write("Ntrans = "+str(ntrans)+"\n")

rmtrx = np.fromfile(f0, dtype = np.int32, count = 9 * ntrans)
f0.read(8)
trans = np.fromfile(f0, dtype = np.double, count = 9 * ntrans)
f0.read(8)
tnp = np.fromfile(f0, dtype = np.double, count = 3 * ntrans)
f0.read(8)
ii=0
if stype[2]==0 :
    M=[]
    for it in range(ntrans) :
        f1.write("M"+str(it)+" =\n")
        formatprint(trans[ii:ii+9], 3, "{:3.0f}")
        M0=[]
        for ix1 in range(3) :
            M1=[]
            for ix2 in range(3) :
                M1.append(trans[ii+ix1*3+ix2])
            M0.append(M1)
        M.append(M0)
        ii+=9
else :
    for it in range(ntrans) :
        f1.write("M"+str(it)+" (r-lattice) =\n")
        formatprint(rmtrx[ii:ii+9], 3, "{:8.4f}")
        f1.write("M"+str(it)+" (real) =\n")
        formatprint(trans[ii:ii+9], 3, "{:3.0f}")
        ii+=9
        f1.write("t"+str(it)+" =\n")
        formatprint(tnp[it*3:it*3+3], 1, "{:8.4f}")
        f1.write()

alatt=np.fromfile(f0, dtype=np.double, count=9)
if stype[2]!=0 :
    f1.write("alatt =\n")
    formatprint(alatt, 3, "{:13.8f}")

invlat=np.fromfile(f0, dtype=np.double, count=9)
f0.read(8)
if stype[2]!=0 :
    f1.write("invlat =\n")
    formatprint(invlat, 3, "{:13.8f}")

chi=np.fromfile(f0, dtype=np.int32, count=ntrans**2)
f0.read(8)
f1.write("chi = \n")
formatprint(chi, ntrans, "{:3}")

chimat=[]
for isym1 in range(ntrans) :
    chimat0=[]
    for isym2 in range(ntrans) :
        chimat0.append(chi[isym1*ntrans+isym2])
    chimat.append(chimat0)

kp_read=np.fromfile(f0, dtype=np.int32, count=3*nwedge)
f0.read(8)
kp=[]  # grid in lattice space
kpr=[] # grid in real space (Bohr)
kprabssum=[]
for irk in range(nwedge) :
    kp.append(kp_read[irk*3:irk*3+3].tolist())
    kpr0=[0.0]*3
    kprabssum0=0.0
    for ix in range(3) :
        kpr0[ix]=(shift[ix]+kp[irk][ix])*gridstep
        kprabssum0+=abs(kpr0[ix])
    kpr.append(kpr0)
    kprabssum.append(kprabssum0)
#(kprabssum_sort,kp_sort,kpr_sort)=zip(*sorted(zip(kprabssum,kp,kpr)))
f1.write("Selected grid points:\n")
f1.write("shift = "+str(shift)+"\n")
f1.write("kpr = (shift(1)+kp)*h\n")
sel_kp=[]
for irk in range(nwedge) :
    #if abs(kpr[irk][0])<1.55 and abs(kpr[irk][1])<1.55 and abs(kpr[irk][2])<1.55  and abs(abs(kpr[irk][0])+abs(kpr[irk][1])+abs(kpr[irk][2])-2.7)<0.05 and abs(abs(kpr[irk][0])-abs(kpr[irk][1]))>0.05:
    if abs(kpr[irk][0])<0.4 and abs(kpr[irk][1])<0.4 and abs(kpr[irk][2])<0.4 :
    #f1.write(f"[{kp_sort[irk][0]:2}{kp_sort[irk][1]:3}{kp_sort[irk][2]:3} ]  [{kpr_sort[irk][0]:5.2f}{kpr_sort[irk][1]:6.2f}{kpr_sort[irk][2]:6.2f} ]\n")
        sel_kp.append(irk)
        f1.write(f"{irk:5d} [{kp[irk][0]:2}{kp[irk][1]:3}{kp[irk][2]:3} ]  [{kpr[irk][0]:5.2f}{kpr[irk][1]:6.2f}{kpr[irk][2]:6.2f} ]\n")
f1.write("\n")

nstate=np.fromfile(f0, dtype=np.int32, count=1)[0]
f0.read(8)
f1.write("nstate = "+str(nstate)+"\n")

irep=np.fromfile(f0, dtype=np.int32, count=nstate)
f0.read(8)
en=np.fromfile(f0, dtype=np.double, count=nstate)
f0.read(8)
occ=np.fromfile(f0, dtype=np.double, count=nstate)
f0.read(8)
en=en*Ry
#irep=np.flip(irep)
for ib in range(nstate) :
    irep[ib]-=1
#en=np.flip(en)
#occ=np.flip(occ)
for ib in range(nstate) :
    if occ[ib]<0.5 :
        nv=ib
        break
nc=nstate-nv
f1.write("   ib irep energy(eV) occ\n")

for ib in range(nv+100,-1,-1) :
#for ib in range(nv+5,nv-7,-1) :
    f1.write(f"{ib+1:6}{irep[ib]:3}{en[ib]:11.6f}{occ[ib]:5.2f}\n")
f1.write("\n")

pot=np.fromfile(f0, dtype=np.double, count=nwedge)
f0.read(8)
chg=np.fromfile(f0, dtype=np.double, count=nwedge)
f0.read(8)
totchg=sum(chg)*gridstep**3*ntrans
f1.write(f"Total charge = {totchg:.2f}\n\n")

f1.write(" === information of k-point 1 ===\n")

# These should be the same as above
nwedge0=np.fromfile(f0, dtype=np.int32, count=1)[0]
ntrans0=np.fromfile(f0, dtype=np.int32, count=1)[0]
f0.read(8)
rmtrx0=np.fromfile(f0, dtype=np.int32, count=9*ntrans0)
f0.read(8)
trans0=np.fromfile(f0, dtype=np.double, count= 9*ntrans0)
f0.read(8)
tnp0=np.fromfile(f0, dtype=np.double, count= 3*ntrans0)
f0.read(8)
alatt0=np.fromfile(f0, dtype=np.double, count=9)
invlat0=np.fromfile(f0, dtype=np.double, count=9)
f0.read(8)
chi0=np.fromfile(f0, dtype=np.int32, count=ntrans0**2)
f0.read(8)
kp_read0=np.fromfile(f0, dtype=np.int32, count=3*nwedge0)
f0.read(8)
nstate0=np.fromfile(f0, dtype=np.int32, count=1)[0]
f0.read(8)
f1.write("nstate0 = "+str(nstate0)+"\n")

indxsave0=np.fromfile(f0, dtype=np.int32, count=nstate0)
f0.read(8)
#f1.write("indxsave0 = "+str(indxsave0)+"\n")

wfn=[]
ncplx=stype[1]+1 # ncplx=2 need to test
for ib in range(nstate0) :
    wfn0=np.fromfile(f0, dtype=np.double, count=nwedge0*ncplx).tolist()
    f0.read(8)
    wfn.append(wfn0)
#f1.write(str(wfn))

# =========================================
# Add any wfn processing code here
# #For example the norm checking: ==========
# f1.write("Example of postprocessing:\n")
# remain=0.0
# for ib in range(nstate0) :
#     summ=0
#     for ig in range(nwedge0) :
#         summ+=wfn[ib][ig]**2
#     remain=max(remain,abs(summ*ntrans0-1))
# f1.write(f"maximum of abs(norm2(wfn[ib])-1): {remain:.2e}\n")
# =========================================

# Checking for symmetries: ================
N=nv+500
deg=[-1]*(N+2)
f1.write("Checking for degeneracy\n")
for ib1 in range(N) :
    if deg[ib1]!=-1 :
        continue
    ideg=0
    for ib2 in range(2,0,-1) :
        if abs(en[ib1]-en[ib1+ib2])<=1e-6 :
            ideg=ib2
            break
    for ib2 in range(ideg+1) :
        deg[ib1+ib2]=ideg

irep_bin={}
irep_new_bin={}

for ib1 in range(N) :
    if irep[ib1] not in irep_bin :
        irep_bin[irep[ib1]]=1
    else :
        irep_bin[irep[ib1]]+=1
    f1.write(f"{ib1+1:5d}{irep[ib1]:2d}{deg[ib1]+1:2d}")
    isigma_d=-1
    for gg in sel_kp :
        if abs(wfn[ib1][gg])>1e-6 :
            isigma_d=1
        f1.write(f"{wfn[ib1][gg]:10.6f}")
    f1.write("\n")
    if deg[ib1]+1==2 :
        irep_new=2 # E
    elif deg[ib1]+1==1 and isigma_d==1 :
        irep_new=0 # A_1
    elif deg[ib1]+1==1 and isigma_d==-1 :
        irep_new=1 # A_2
    elif deg[ib1]+1==3 and isigma_d==-1 :
        irep_new=3 # T_1
    elif deg[ib1]+1==3 and isigma_d==1 :
        irep_new=4 # T_2
    if irep_new not in irep_new_bin :
        irep_new_bin[irep_new]=1
    else :
        irep_new_bin[irep_new]+=1

f1.write("Old irep count: "+str(irep_bin)+"\n")
f1.write("New irep count: "+str(irep_new_bin)+"\n")

#f1.write("Checking for symmetries <i|j>, only reduced BZ\n")
#for ib1 in range(N) :
#    for ib2 in range(N) :
#        summ=0
#        for ig in range(nwedge0) :
#            summ+=wfn[ib1][ig]*wfn[ib2][ig]
#        if abs(summ)<=1e-6 : 
#            summ=0
#        f1.write(str(ib1)+" "+str(ib2)+" "+str(summ*ntrans0)+"\n")
#
## creating a map of all kpr in full BZ
#f1.write("Checking for symmetries <i|j>, FULL BZ\n")
#
#kpr_full=[]
#kmap=[]
#kfactor=[]
#for isym1 in range(ntrans) :
#    kpr0=[]
#    kmap0=[]
#    kfactor0=[]
#    for isym2 in range(ntrans) :
#        for ig1 in range(nwedge0) :
#            # calculate the new kpr in full BZ, test if duplicate, then write
#            kpr1=v3tm3(kpr[ig1],M[isym2])
#            #fdup=False
#            #for ig2 in range(len(kpr0)) :
#            #    summ=0
#            #    for ix in range(3) :
#            #        summ+=abs(kpr0[ig2][ix]-kpr1[ix])
#            #    if summ<1e-8 :
#            #        fdup=True
#            #if fdup :
#            #    continue
#            kpr0.append(kpr1)
#            kmap0.append(ig1)
#            kfactor0.append(chimat[isym1][isym2])
#
#    kpr_full.append(kpr0)
#    kmap.append(kmap0)
#    kfactor.append(kfactor0)
#
#f1.write("kpr: \n")
#f1.write(str(len(kpr))+'\n')
#f1.write("kpr_full: \n")
#f1.write(str(len(kpr_full[0]))+'\n')
##f1.write(str(kmap)+'\n')
##f1.write(str(kfactor)+'\n')
#
#nwdge0_full=len(kpr_full[0])
#f1.write("nwdge0_full = "+str(nwdge0_full)+"\n")
#N=10
#f1.write("info of first "+str(N)+" bands:\n")
#for ib in range(N) :
#    f1.write(str(ib)+" "+str(irep[ib])+"\n")
#for ib1 in range(N) :
#    for ib2 in range(N) :
#        summ=0
#        for ig in range(nwdge0_full) :
#            ir1=irep[ib1]
#            ir2=irep[ib2]
#            ig1=kmap[irep[ib1]][ig]
#            ig2=kmap[irep[ib2]][ig]
#
#            #coul=0
#            ## coul: 1/|r-r'|
#            #for ix in range(3) :
#            #    coul=(kpr_full[ib1][ig1][ix]-kpr_full[ib2][ig2][ix])**2
#            #coul=coul**(-0.5)
#            summ+=wfn[ib1][ig1]*wfn[ib2][ig2]*kfactor[ir1][ig]*kfactor[ir2][ig]
#        #if abs(summ)<=1e-8 :
#        #    summ=0
#        f1.write(str(ib1)+" "+str(irep[ib1])+" "+str(ib2)+" "+str(irep[ib2])+" "+str(summ)+"\n")

f0.close()
f1.close()
os.system("cat wfninfo.out")
