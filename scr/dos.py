#!/bin/python3

# print DOS

# Input: INCAR, DOSCAR

import sys
import matplotlib.pyplot as plt

f0=open("INCAR","r")
line=f0.readlines()
f0.close()

ispin=1
nedos=301
for l in range(len(line)) :
    word=line[l].replace("="," ").split()
    if len(word)==0 or word[0][0]=="!" or word[0][0]=="#" :
        continue
    if word[0]=="ISPIN":
        ispin=int(word[1])
    if word[0]=="NEDOS":
        nedos=int(word[1])

if ispin==1:
    print("only support spin=2")
    sys.exit()

f1=open("DOSCAR","r")
line=f1.readlines()
f1.close()

x=[]
dos=[[],[]]
ef=float(line[5].split()[2])
for l in range(6,6+nedos) :
    word=line[l].split()
    x.append(float(word[0]))
    dos[0].append(float(word[1]))
    dos[1].append(-float(word[2]))

fig, ax = plt.subplots()
ax.axvline(x=ef,color="black",linewidth=0.75)
ax.plot(x,dos[0],color="tab:blue",linewidth=1)
ax.plot(x,dos[1],color="tab:orange",linewidth=1)
ax.set_xlim([0,15])
ax.set_ylim([-200,200])
ax.set_ylabel("DOS")
ax.set_xlabel("Energy (eV)")
fig.savefig("dos.png",dpi=300,bbox_inches='tight')
plt.show()
