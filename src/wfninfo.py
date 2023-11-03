#!/usr/bin/env python3
# This tool prints some information from parsec wavefunction file parsec.dat
# Now only support the real wavefunction of cluster, which has 1 k-point (Gamma)

# python3 wfninfo.py

# outputs (wfninfo.out) are shown on the screen

import sys
import os
import numpy as np

def formatprint(mat, dim0, fstr) :
    # dim0 is the number of rows
    # example of fstr: "{:3}" or "{:.2f}"
    dim1=len(mat)//dim0
    strformat=dim1*fstr
    for ii in range(dim0) :
        f1.write(strformat.format(*mat[ii*dim1:(ii+1)*dim1]))
        f1.write("\n")

Ry=13.605703976

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

ndim=np.fromfile(f0, dtype=np.int32, count=1)[0]
f0.read(8)
f1.write("Ndim = "+str(ndim)+"\n")

shift=np.fromfile(f0, dtype=np.double, count=3)
f0.read(8)

nwedge=np.fromfile(f0, dtype=np.int32, count=1)[0]
ntrans=np.fromfile(f0, dtype=np.int32, count=1)[0]
f0.read(8)
f1.write("Nwedge = "+str(nwedge)+" (Ndim in irreducible wedge)\n")
f1.write("Ntrans = "+str(ntrans)+"\n")

rmtrx=np.fromfile(f0, dtype=np.int32, count=9*ntrans)
f0.read(8)
trans=np.fromfile(f0, dtype=np.double, count= 9*ntrans)
f0.read(8)
tnp=np.fromfile(f0, dtype=np.double, count= 3*ntrans)
f0.read(8)
ii=0
if stype[2]==0 :
    for it in range(ntrans) :
        f1.write("M"+str(it)+" =\n")
        formatprint(trans[ii:ii+9], 3, "{:3.0f}")
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
(kprabssum_sort,kp_sort,kpr_sort)=zip(*sorted(zip(kprabssum,kp,kpr)))
f1.write("kp of top 10 grid points in wedge\n")
f1.write("shift = "+str(shift)+"\n")
f1.write("kpr = (shift(1)+kp)*h\n")
for irk in range(10) :
    f1.write(f"[{kp_sort[irk][0]:2}{kp_sort[irk][1]:3}{kp_sort[irk][2]:3} ]  [{kpr_sort[irk][0]:5.2f}{kpr_sort[irk][1]:6.2f}{kpr_sort[irk][2]:6.2f} ]\n")
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
irep=np.flip(irep)
en=np.flip(en)
occ=np.flip(occ)
for ib in range(nstate) :
    if occ[ib]>0.5 :
        nc=ib
        break
nv=nstate-nc
f1.write("   ib irep energy(eV) occ\n")
for ib in range(nc-6,nc+6) :
    f1.write(f"{nstate-ib:6}{irep[ib]:3}{en[ib]:11.6f}{occ[ib]:5.2f}\n")
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

# Add any wfn processing code here
# For example the norm checking:
f1.write("Example of postprocessing:\n")
remain=0.0
for ib in range(nstate0) :
    summ=0
    for ig in range(nwedge0) :
        summ+=wfn[ib][ig]**2
    remain=max(remain,abs(summ*ntrans0-1))
f1.write(f"maximum of abs(norm2(wfn[ib])-1): {remain:.2e}\n")

f0.close()
f1.close()
os.system("cat wfninfo.out")
