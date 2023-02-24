#!/usr/bin/env python3

# generate the k grid for wannier90 and corresponding qe pw.x nscf.

# python3 wankgrid.py Nx Ny Nz

# Output: wannier90_kgrid.dat qe_kgrid.dat

import sys

if len(sys.argv)<4 :
    print("Need more arguments:")
    print("python3 wankgrid.py Nx Ny Nz")
    sys.exit()

Nx=int(sys.argv[1])
Ny=int(sys.argv[2])
Nz=int(sys.argv[3])
N=Nx*Ny*Nz
w=1.0/N

x=[0.0]*Nx
y=[0.0]*Ny
z=[0.0]*Nz

for ix in range(Nx) :
    x[ix]=1.0*ix/Nx
for iy in range(Ny) :
    y[iy]=1.0*iy/Ny
for iz in range(Nz) :
    z[iz]=1.0*iz/Nz

f1=open("wannier90_kgrid.dat","w")
f2=open("qe_kgrid.dat","w")

f1.write("mp_grid =  "+str(Nx)+"  "+str(Ny)+"  "+str(Nz)+"\n\n")
f1.write("begin kpoints\n")
f2.write("K_POINTS crystal\n")
f2.write(str(N)+"\n")
for iz in range(Nz) :
    for iy in range(Ny) :
        for ix in range(Nx) :
            f1.write("  {:.8f}".format(x[ix]))
            f1.write("  {:.8f}".format(y[iy]))
            f1.write("  {:.8f}\n".format(z[iz]))
            f2.write("  {:.8f}".format(x[ix]))
            f2.write("  {:.8f}".format(y[iy]))
            f2.write("  {:.8f}".format(z[iz]))
            f2.write("  {:.8f}\n".format(w))
f1.write("end kpoints\n\n")
f2.write("\n")

f1.close()
f2.close()
