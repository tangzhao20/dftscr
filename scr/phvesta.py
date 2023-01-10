#!/bin/env python3

# Read the vibration vectors and add to the vesta files

# Input: qpoints.yaml, structure.vesta

import yaml

scal=4
adia=0.2  # diameter of the arrow
rgb=[0,128,255] # arrow colors



with open("qpoints.yaml") as f0:
    qp=yaml.full_load(f0)

with open("structure.vesta") as f0:
    line=f0.readlines()



Nl=len(line)
Nb=len(qp["phonon"][0]["band"])
Na=len(qp["phonon"][0]["band"][0]["eigenvector"])

freq=[]
vec=[]
for ib in range(Nb):
    freq.append(qp["phonon"][0]["band"][0]["frequency"])
    vec0=[]
    for ia in range(Na) :
        vec1=[0]*3
        for id in range(3) :
            vec1[id]=qp["phonon"][0]["band"][ib]["eigenvector"][ia][id][0]*scal
        vec0.append(vec1)
    vec.append(vec0)
    

# writting from here

for ib in range(Nb) :
    il=0
    f1=open("ph_"+str(ib+1)+".vesta","w")
    while il<Nl:
        # copy the original vesta files
        if line[il]=="VECTR\n" :
            break
        f1.write(line[il])
        il+=1
    # write the vectors
    f1.write("VECTR\n")
    for ia in range(Na):
        f1.write(str(ia+1)+" "+str(vec[ib][ia][0])+" "+str(vec[ib][ia][1])+" "+str(vec[ib][ia][2])+"\n")
        f1.write(str(ia+1)+" 0 0 0 0\n")
        f1.write("0 0 0 0 0\n")
    f1.write("0 0 0 0 0\n")
    f1.write("VECTT\n")
    for ia in range(Na):
        f1.write(str(ia+1)+" "+str(adia)+" "+str(rgb[0])+" "+str(rgb[1])+" "+str(rgb[2])+" 0\n")
    f1.write("0 0 0 0 0\n")
    il=il+4
    while il<Nl:
        # copy the original vesta files
        f1.write(line[il])
        il+=1
    f1.close()


#print(vec["phonon"][0]["band"][0]["eigenvector"])
