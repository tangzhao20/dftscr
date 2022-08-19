eigenval.py :

eigenval.py (E1) (E2)
Making the band structure plot

if E2 exists, the energy range is (E1,E2)
if only E1 exists, the energy range is (-E1,E1)
if neither exists, the energy range is (-5 eV, 5 eV)

Input: EIGENVAL, KPOINTS, POSCAR (DOSCAR)

########################################

proj.py :

proj.py atoms orbitals (E1) (E2)
Make the projected band structure plot

if E2 exist, the energy range is (E1,E2)
if only E1 exists, the energy range is (-E1,E1)
if neither exists, the energy range is (-5 eV, 5 eV)

Input: POSCAR, EIGENVAL, KPOINTS, PROCAR, (DOSCAR)

########################################

dos.py :

print DOS

Input: INCAR, DOSCAR

########################################

toten_fit.py :

print the toten vs. distortion coordinate 

Input: POSCAR_i, pos_*/OSZICAR, pos_*/POSCAR, (pos_*/EIGENVAL)

########################################

phproj.py :

Calculate the projection of a distortin onto the phonon modes,
and use it as the weight to average the frequency

Input: qpoints.yaml, POSCAR_i, POSCAR_f

########################################

ZT
