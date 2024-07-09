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
python posconvert.py package1 package2 (filename1)
```
Support packages: VASP, QE (.in or .xml), Paratec, Parsec, .xyz

VASP input: POSCAR  
QE input: nscf.in/scf.in/relax.in  
QExml input: \*.xml  
Paratec input: input  
PARSEC input: parsec.in  
.xyz input: \*.xyz  

Optional input: posconvert.in  
If the posconvert.in file exist, this code modify to the structure before output. The supported opreation includes move, rotate, flip, vacuum, supercell, add atom, or delete atom. An example of [posconvert.in](data/inputs/posconvert.in) can be find in [data/inputs](data/inputs) directory.

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
python afm.py (vasp)
```
Input: afm.in, tip.xyz, sample.parsec\_st.dat, (\*\_POTRE.DAT)  
Output : parsec\_st\_\*\_\*.dat, (parsec\_st\_spot.dat,) manual\_\*\_\*.dat, steps.dat  
An example of [afm.in](data/inputs/afm.in) can be find in [data/inputs](data/inputs) directory. If fdet is set in afm.in, parsec\_st\_spot.dat will be generated, and the pseudopotential files (\*\_POTRE.DAT) will be read to calculate Nb. Add `vasp` at the end to write an example structure (tip + sample) in VASP format.

The Bash script can be used to create directories and prepare files by
```
afm.sh seq 
```
Input: afm.in, job.sh, parsec.in.head, parsec\_st\_\*\_\*.dat, (parsec\_st\_spot.dat,) manual\_\*\_\*.dat, steps.dat  
Output: seq\_\*\_\*, (spot,) \*/parsec.in, \*/job.sh, sbatch.log  

If the fdet option is set, enter the spot directory and submit the job manually before iteratively submitting all other jobs. If fdet is not set, use the following command to submit the jobs:
```
afm.sh sbatch
```

After all the calculations are done, use the python script [afmplot.py](src/afmplot.py) to make the plots.  
```
python afmplot.py (iz) (atom) (tilt) (bohr)
```
Input: afm.in, steps.dat, (toten.dat or seq\_\*\_\*/parsec.out)  
Output: afm\_\*.png, (toten.dat)  
In the first round, this code read the total energies from seq\_\*\_\*/parsec.out files and write to toten.dat. After that, the toten.dat will be read.  
`iz` represents the index of layers to be calculated. For the simple scenario of computing 3 z values for the tip, `iz` should be set to 2 as default. Add the `atom` option to display the atom positions. Only the top layer within 1 A are plotted. Add the `tilt` option to use the tilt correction. The tilt correction method is on development for now. Add the `bohr` option to use the Bohr as the length unit.  

---

ZT
