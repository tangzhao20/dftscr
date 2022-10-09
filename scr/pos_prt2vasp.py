import sys

f1=open("input","r")
line=f1.readlines()
f1.close()

Flc=False
Fat=False
lc=[]
attype=[]
atom=[]
iatom=[]
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
            lc.append([float(word[1]),float(word[2]),float(word[3])])
        elif len(word)>=2 and word[0]=="volume" :
            volume=float(word[1])
    elif Fat==True :
        if len(word)>=2 and word[0]=="newtype" :
            attype.append(word[1])
            iatom.append(0)
        elif len(word)>=4 and word[0]=="coord" :
            atom.append([float(word[1]),float(word[2]),float(word[3])])
            iatom[len(iatom)-1]=iatom[len(iatom)-1]+1

detlc=(lc[0][0]*(lc[1][1]*lc[2][2]-lc[1][2]*lc[2][1])+lc[0][1]*(lc[1][2]*lc[2][0]-lc[1][0]*lc[2][2])+lc[0][2]*(lc[1][0]*lc[2][1]-lc[1][1]*lc[2][0]))
factor=(volume/6.748334503468005/detlc)**(1.0/3.0)

for i in range(3) :
    for j in range(3) :
        lc[i][j]=lc[i][j]*factor

for i in range(3) :
    for j in range(3) :
        if(lc[i][j])==-0.0 :
            lc[i][j]=0.0
for i in range(len(atom)) :
    for j in range(3) :
        if(atom[i][j])==-0.0 :
            atom[i][j]=0.0

f2=open("POSCAR","w")
if len(sys.argv)==1 :
    f2.write("SYSTEM\n1.0\n")
else :
    f2.write(sys.argv[1]+"\n1.0\n")
for i in range(3) :
    f2.write("  "+str(lc[i][0])+"  "+str(lc[i][1])+"  "+str(lc[i][2])+"\n")
for i in range(len(attype)) :
    f2.write("   "+attype[i])
f2.write("\n")
for i in range(len(iatom)) :
    f2.write("   "+str(iatom[i]))
f2.write("\nDirect\n")
for i in range(len(atom)) :
    f2.write("  "+str(atom[i][0])+"  "+str(atom[i][1])+"  "+str(atom[i][2])+"\n")
f2.close()
