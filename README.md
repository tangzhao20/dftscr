## `bands.py`

To plotting the band structure plot:
```
python bands.py package (E1) (E2)
```
if `E2` exists, the energy range is [`E1`,`E2`]  
if only `E1` exists, the energy range is [-`E1`,`E1`]  
if neither exists, the energy range is [-5 eV, 5 eV]

Support `package`s: `vasp`, `qe`  
VASP inputs: `EIGENVAL`, `KPOINTS`, `POSCAR`, (`DOSCAR`)  
QE inputs: `*.xml`, `nscf.in`, `kpath.in`

To plot projected band structure from VASP:
```
python bands.py vaspproj atoms orbitals (E1) (E2)
```
Inputs: `EIGENVAL`, `KPOINTS`, `POSCAR`, `PROCAR`, (`DOSCAR`)

To plot band structure from wannier90:
```
python bands.py wan (pad=*) (E1) (E2)
```
Inputs: `nscf.in`, `../bands/*.xml`, `*_band.kpt`, `*_band.dat`, `kpath.in`

---

## `posconvert.py`

To convert the format of the crystal structure.
```
python posconvert.py package1 package2
```
Support `package`s: `vasp`, `qe`, `paratec`, `parsec`

VASP input: `POSCAR`  
QE input: `scf.in`  
Paratec input: `input`

---

## `kconvert.py`

To convert the format of the k-point path.
```
python kconvert.py package1 package2 (N)
```
Support `package`s: `vasp`, `qe`, `kpathin`

VASP input: `KPOINTS`  

Create nscf k-point path for QE band structure calculations
```
python kconvert.py kpathin qe (N)
```

---

## `dos.py`

To print the DOS  
Input: `INCAR`, `DOSCAR`

---

## `toten_fit.py`

To print the total energy vs. distortion coordinate  
Input: `POSCAR_i`, `pos_*/OSZICAR`, `pos_*/POSCAR`, (`pos_*/EIGENVAL`)

---

## `phproj.py`

To calculate the projection of a distortin onto the phonon modes, and use it as weight to average the frequency  
Input: `qpoints.yaml`, `POSCAR_i`, `POSCAR_f`

---

## `wavecar.py`

To read the KS orbitals from the `WAVECAR`, and write into seperate files for VESTA plot:
```
python wavecar.py ik ib ispin
```
Input: `kpath.in`

Examples of kpath.in can be find in `data` directory

---

ZT
