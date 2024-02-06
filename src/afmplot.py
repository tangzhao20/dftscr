#!/usr/bin/env python3

# This script reads the afm simulation results and makes the plot

import sys
import os
from commons import load_palette
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ==================== read the input file ====================
x_spacing=0.6
y_spacing=0.6
z_sampling=[5.7,6.0,6.3]
z_spacing=0.3
f1=open("afm.in","r")
line=f1.readlines()
f1.close()
for l in line :
    word=l.split()
    if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
        continue
    if word[0]=="x_range" :
        x_range=[float(word[1]),float(word[2])]
    elif word[0]=="y_range" :
        y_range=[float(word[1]),float(word[2])]
    elif word[0]=="x_spacing" :
        x_spacing=float(word[1])
    elif word[0]=="y_spacing" :
        y_spacing=float(word[1])
    elif word[0]=="z_sampling" :
        z_sampling=word[1:]
        for iz in range(range(z_sampling)) :
            try :
                z_sampling[iz]=float(z_sampling[iz])
            except :
                del zsampling[iz:]
                break
    elif word[0]=="parallel" :
        parallel=int(word[1])
    else :
        print("Warning: keyword "+word[0]+" is not defined.")

x=x_range[0]
nx=0
xlist=[]
while (x<x_range[1]+1e-6) :
    xlist.append(x)
    x+=x_spacing
    nx+=1
y=y_range[0]
ny=0
ylist=[]
while (y<y_range[1]+1e-6) :
    ylist.append(y)
    y+=y_spacing
    ny+=1

f3=open("steps.dat","r")
line=f3.readlines()
f3.close()
steplist=[]
for l in line :
    word=l.split()
    if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
        continue
    steplist.append(int(word[0]))

lxincrease=True
k=0
ip=0
# movelist here is the index (ix,iy) of each point
movelist=[]
for iy in range(ny) :
    if lxincrease :
        for ix in range(nx) :
            if k==0 :
                movelist.append([])
            movelist[-1].append([ix,iy])
            k+=1
            if k==steplist[ip] :
                ip+=1
                k=0
    else :
        for ix in range(nx-1,-1,-1) :
            if k==0 :
                movelist.append([])
            movelist[-1].append([ix,iy])
            k+=1
            if k==steplist[ip] :
                ip+=1
                k=0
    lxincrease=not lxincrease



files = os.listdir()
fkts=False
if "kts.dat" in files :
    f4=open("kts.dat","r")
    line=f4.readlines()
    f4.close()
    il=0
    kts=[]
    for ix in range(nx) :
        kts0=[]
        for iy in range(ny) :
            kts0.append(float(line[il].split()[2]))
            il+=1
        kts.append(kts0)

else:
    # initialize toten matrix
    toten=[]
    for iz in range(3) :
        toten0=[]
        for iy in range(ny) :
            toten0.append([0.0]*nx)
        toten.append(toten0)

    for iz in range(3) :
        for ip in range(parallel) :
            istep=-1
            filename="seq_"+str(iz+1)+"_"+str(ip+1)+"/parsec.out"
            f1=open(filename,"r")
            line=f1.readlines()
            f1.close()
            
            for l in line: 
                word=l.split()
                if len(word)==0 or word[0][0]=="#" or word[0][0]=="!" :
                    continue
                if len(word)>=2 and word[0]=="Starting" and word[1]=="SCF..." :
                    istep+=1
                if len(word)>=5 and word[0]=="Total" and word[1]=="Energy" and word[2]=="=" :
                    toten[iz][movelist[ip][istep][1]][movelist[ip][istep][0]]=float(word[3]) # in Ry


    kts=[]
    for ix in range(nx) :
        kts0=[]
        for iy in range(ny) :
            kts1=(0.25*toten[0][ix][iy]-0.5*toten[1][ix][iy]+0.25*toten[2][ix][iy])/z_spacing**2
            kts0.append(kts1)
        kts.append(kts0)

    f4=open("kts.dat","w")
    for ix in range(nx) :
        for iy in range(ny) :
            f4.write(str(ix+1)+" "+str(iy+1)+" "+str(kts[ix][iy])+"\n")
    f4.close()

palette=load_palette() #color palette: [blue, orange, gray, white, black]
mpl.rcParams["font.sans-serif"].insert(0,"Noto Sans")
mpl.rcParams.update({'font.size': 14})

plt.imshow(kts,interpolation='gaussian',cmap=cm.gray,extent=[x_range[0],x_range[1],y_range[0],y_range[1]])

plt.colorbar()

plt.savefig("afm.png",dpi=600)
