from vasp import filein 
import sys
import os
import math

theta=math.radians(float(sys.argv[1]))
costheta=math.cos(theta)
sintheta=math.sin(theta)

poscar1=filein.POSCAR()
aa=poscar1.lc[0][0]
bb=poscar1.lc[0][1]
cc=poscar1.lc[1][0]
dd=poscar1.lc[1][1]
poscar1.lc[0][0]=aa*costheta-bb*sintheta
poscar1.lc[0][1]=aa*sintheta+bb*costheta
poscar1.lc[1][0]=cc*costheta-dd*sintheta
poscar1.lc[1][1]=cc*sintheta+dd*costheta

poscar1.filewrite("POSCAR.new")
