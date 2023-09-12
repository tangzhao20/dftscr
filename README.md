## [bands.py](src/bands.py)

To plotting the band structure plot:
```
python bands.py package (E1) (E2)
```
if E2 exists, the energy range is [E1, E2]  
if only E1 exists, the energy range is [-E1, E1]  
if neither exists, the energy range is [-5 eV, 5 eV]

Support packages: VASP, QE, Parsec  
VASP inputs: EIGENVAL, KPOINTS, POSCAR, (DOSCAR)  
QE inputs: \*.xml, kpath.in  
Parsec inputs: bands.dat, parsec.in, kpath.in

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

Examples of kpath.in can be find in [data/kpath](data/kpath) directory

---

## [silicon.py](src/silicon.py)

This code creates the molecular structure of silicon clusters, given the radius in angstrom.
```
python silicon.py radius
```

---

## [afm.py](src/afm.py)

This Python code, along with the Bash script [afm.sh](src/afm.sh), prepares the job directories and files for the AFM simulation.

The Python code can be used to prepare the structure files by
```
python afm.py
```
Input: afm.in, tip.xyz, sample.xyz  
Output : parsec\_st.dat, manual\_\*\_\*.dat  
The structure files manual\_\*\_\*.dat will be moved into the seq\_\*\_\* directories in the Bash script.

The Bash script can be used to create directories and prepare files by
```
dftscr_dir/jobs/afm.sh
```
Input: afm.in, tip.xyz, sample.xyz, job.sh, parsec.in.head  
Output: seq\_\*\_\*, seq\_\*\_\*/parsec.in, seq\_\*\_\*/job.sh

Examples of [afm.in](data/inputs/afm.in) can be find in [data/inputs](data/inputs) directory.

---

## [posvacuum.py](src/posvacuum.py)

This script adds the vacuum layer to the structure.
```
python posvacuum.py z_vac
```
Input: POSCAR  
Output: POSCAR.new

---

ZT
