## [bands.py](src/bands.py)

To plotting the band structure plot:
```
python bands.py package (E1) (E2)
```
if E2 exists, the energy range is [E1, E2]  
if only E1 exists, the energy range is [-E1, E1]  
if neither exists, the energy range is [-5 eV, 5 eV]

Support packages: VASP, QE  
VASP inputs: EIGENVAL, KPOINTS, POSCAR, (DOSCAR)  
QE inputs: \*.xml, kpath.in

To plot projected band structure:
```
python bands.py package atoms orbitals (E1) (E2)
```
Support packages: VASP (vaspproj), QE (qeproj)  
VASP inputs: EIGENVAL, KPOINTS, POSCAR, PROCAR, (DOSCAR)  
QE inputs: \*.xml, kpath.in, projwfc.out

To plot band structure from wannier90:
```
python bands.py wan (pad=*) (E1) (E2)
```
Inputs: nscf.in, ../bands/\*.xml, \*\_band.kpt, \*\_band.dat, kpath.in

---

## [posconvert.py](src/posconvert.py)

To convert the format of the crystal structure.
```
python posconvert.py package1 package2
```
Support packages: VASP, QE (.in or .xml), Paratec, Parsec, .xyz

VASP input: POSCAR
QE input: scf.in
Paratec input: input

---

## [kconvert.py](src/kconvert.py)

To convert the format of the k-point path.
```
python kconvert.py package1 package2 (N)
```
Support packages: VASP, QE, kpathin

VASP input: KPOINTS  

Create nscf k-point path for QE band structure calculations
```
python kconvert.py kpathin qe (N)
```

---

## [dos.py](src/dos.py)

To plot DOS  
```
python dos.py (v) package (E1) (E2)
```
Support packages: VASP, QE

VASP input: DOSCAR  
QE input: \*.dos

To plot projected DOS:
```
python dos.py (v) package atoms orbitals (E1) (E2)
```
Support packages: QE (qeproj)  
QE input: \*.dos \*.xml \*.pdos\_atm#\*(\*)\_wfc#\*(\*)

---

## [toten\_fit.py](src/toten_fit.py)

To print the total energy vs. distortion coordinate  
Input: POSCAR\_i, pos\_\*/OSZICAR, pos\_\*/POSCAR, (pos\_\*/EIGENVAL)

---

## [phproj.py](src/phproj.py)

To calculate the projection of a distortin onto the phonon modes, and use it as weight to average the frequency  
Input: qpoints.yaml, POSCAR\_i, POSCAR\_f

---

## [wavecar.py](src/wavecar.py)

To read the KS orbitals from the WAVECAR, and write into seperate files for VESTA plot:
```
python wavecar.py ik ib ispin
```
Input: kpath.in

Examples of kpath.in can be find in [data](data/kpath) directory

---

## [silicon.py](src/silicon.py)

This code creates the silicon, given the radius in angstrom.
```
python silicon.py radius
```

---
ZT
