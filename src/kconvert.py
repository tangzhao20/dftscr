#!/usr/bin/env python3

# create nscf k-point path for QE band structure calculations
# python kconvert.py package1 package2 (N)

# input: kpath.in

import sys
from classes import kpoints_band
from commons import load_packagename

packagename=load_packagename()

if len(sys.argv)<3:
    print("Need more arguments.")
    print("python kconvert.py package1 package2 (N)")
    sys.exit()
elif len(sys.argv)==3 :
    N=40
else :
    N=int(sys.argv[1])

package1=sys.argv[1]
package2=sys.argv[2]

if package1 in packagename["vasp"] :
    kpoints1=kpoints_band.KPOINTS_band()
elif package1 in packagename["kpathin"] :
    kpoints1=kpoints_band.KPOINTS_band(empty=True)
    kpoints1.fileread_kpathin(nk_line=N)
else :
    print("Package "+package1+" input is not supported yet.")
    print("python3 kconvert.py package1 package2 (N)")
    sys.exit()


if package2 in packagename["kpathin"] :
    kpoints1.filewrite_kpathin()
elif package2 in packagename["qe"] :
    kpoints1.filewrite_qe()
elif package2 in packagename["vasp"] :
    kpoints1.filewrite_vasp()
else :
    print("Package "+package2+" input is not supported yet.")
    print("python3 kconvert.py package1 package2 (N)")
    sys.exit()
