# DFT scripts

This collection of tools is designed to analyze and visualize data from first-principles electronic structure calculations. It includes features such as plotting band structures, converting atomic structures, analyzing the density of states (DOS), and preparing inputs for calculations. Various commonly used calculation packages are supported.

## Installation
Clone the repository and use it as is. The Python scripts are written in Python 3. For convenience, you may add the [src](src) directory to your PATH. To install the required dependencies, use the provided [requirements.txt](requirements.txt) file:
```bash
pip install -r requirements.txt
```

## Band structure plotting
**Script:** [bands.py](src/bands.py)  

For basic band structure plotting, use:  
```bash
python bands.py package (E1) (E2)
```
* If both `E1` and `E2` are provided, the energy range is [E1, E2] (in eV).  
* If only `E1` is provided, the energy range is [-E1, E1] (in eV).  
* If neither is provided, the energy range is [-5 eV, 5 eV].  

**Supported formats:** VASP, QE, PARSEC  

**VASP inputs:** EIGENVAL, KPOINTS, POSCAR, (DOSCAR)  
**QE inputs:** \*.xml, kpath.in  
**PARSEC inputs:** bands.dat, parsec.in, kpath.in

To plot projected band structures:  
```bash
python bands.py package atoms orbitals (E1) (E2)
```
**Supported formats:** VASP (vaspproj), QE (qeproj)  

**VASP inputs:** EIGENVAL, KPOINTS, POSCAR, PROCAR, (DOSCAR)  
**QE inputs:** \*.xml, kpath.in, projwfc.out

To plot band structures from Wannier90:  
```bash
python bands.py wan (pad=*) (E1) (E2)
```
**Inputs:** nscf.in, ../bands/\*.xml, \*\_band.kpt, \*\_band.dat, kpath.in

## Atomic structure conversion
**Script:** [posconvert.py](src/posconvert.py)  

To convert the format of the atomic structure:  
```bash
python posconvert.py package1 package2 (filename1)
```
**Supported formats:** VASP, QE (.in or .xml), Paratec, PARSEC, .xyz  

**VASP input:** POSCAR  
**QE input:** nscf.in/scf.in/relax.in  
**QExml input:** \*.xml  
**Paratec input:** input  
**PARSEC input:** parsec.in  
**.xyz input:** \*.xyz  

**Optional input:** posconvert.in  

If the optional input file posconvert.in exists, this code modifies the structure before output. The supported operations include move, rotate, flip, vacuum, supercell, add atom, or delete atom. An example of [posconvert.in](data/inputs/posconvert.in) can be found in [data/inputs](data/inputs) directory.

## *k*-point path conversion  
**Scripts:** [kconvert.py](src/kconvert.py)

To convert the format of the *k*-point path:  
```bash
python kconvert.py package1 package2 (N)
```
**Supported formats:** VASP, QE, kpath.in

**VASP input:** KPOINTS  

Examples of the kpath.in file can be found in the [data/kpaths](data/kpaths) directory.

A common use of this code is to create the *k*-point path for QE band structure calculations:  
```bash
python kconvert.py kpathin qe (N)
```

## DOS plotting
**Script:** [dos.py](src/dos.py)  

To plot DOS:  
```bash
python dos.py (v) package (E1) (E2)
```
* `v` indicates plot vertically   
* `E1` and `E2` define the energy range, as described in the band structure section

**Supported formats:** VASP, QE  

**VASP input:** DOSCAR  
**QE input:** \*.dos

To plot projected DOS:  
```bash
python dos.py (v) package atoms orbitals (E1) (E2)
```
**Support format:** QE (qeproj)  

**QE inputs:** \*.dos \*.xml \*.pdos\_atm#\*(\*)\_wfc#\*(\*)

## Total energy fitting
**Script:** [toten\_fit.py](src/toten_fit.py)  

This script prints the total energy vs. distortion coordinate.  
**Inputs:** POSCAR\_i, pos\_\*/OSZICAR, pos\_\*/POSCAR, (pos\_\*/EIGENVAL)

## Phonon projection
**Script:** [phproj.py](src/phproj.py)

This script calculates distortion projection onto phonon modes and uses it as the weight to average the frequency.  
**Inputs:** qpoints.yaml, POSCAR\_i, POSCAR\_f

## KS orbitals extraction
**Script:** [wavecar.py](src/wavecar.py)

This script reads the KS orbitals from the WAVECAR, and writes into separate files for the VESTA plot:  
```bash
python wavecar.py ik ib ispin
```
**Input:** WAVECAR, POSCAR

## Silicon nanocrystals creation
**Scripts:** [silicon.py](src/silicon.py), [sidef.py](src/sidef.py)

[silicon.py](src/silicon.py) creates the atomic structure of silicon nanocrystals, given the radius in angstrom:  
```bash
python silicon.py radius
```

[sidef.py](src/sidef.py) prints information of atoms in a silicon nanocrystal:
```bash
python3 sidef.py filename
```

## AFM simulation preparation
**Scripts:** [afm.py](src/afm.py), [afm.sh](src/afm.sh), [afmplot.py](src/afmplot.py)

[afm.py](src/afm.py) and [afm.sh](src/afm.sh) prepare the files and job directories for the AFM simulation. The Python code can be used to prepare the structure files by
```bash
python afm.py (vasp)
```

**Inputs:** afm.in, tip.xyz, sample.parsec\_st.dat, (\*\_POTRE.DAT)  

**Outputs:** parsec\_st\_\*\_\*.dat, (parsec\_st\_spot.dat,) manual\_\*\_\*.dat, steps.dat  

An example of [afm.in](data/inputs/afm.in) can be found in the [data/inputs](data/inputs) directory. If `fdet` is set in afm.in, parsec\_st\_spot.dat will be generated, and the pseudopotential files \*\_POTRE.DAT will be read to calculate Nb. Add `vasp` at the end to write an example structure (tip + sample) in VASP format.

The Bash script can be used to create directories and prepare files by
```bash
afm.sh seq 
```

**Inputs:** afm.in, job.sh, parsec.in.head, parsec\_st\_\*\_\*.dat, (parsec\_st\_spot.dat,) manual\_\*\_\*.dat, steps.dat  

**Outputs:** seq\_\*\_\*, (spot,) \*/parsec.in, \*/job.sh, sbatch.log  

If the `fdet` option is set, enter the spot directory and submit the job manually before iteratively submitting all other jobs. If `fdet` is not set, use the following command to submit the jobs:
```bash
afm.sh sbatch
```

After all the calculations are done, use the Python script [afmplot.py](src/afmplot.py) to make the plots:  
```bash
python afmplot.py (iz) (atom) (tilt) (bohr)
```

**Inputs:** afm.in, steps.dat, (toten.dat or seq\_\*\_\*/parsec.out)  

**Outputs:** afm\_\*.png, (toten.dat)  

In the first round, this code reads the total energies from seq\_\*\_\*/parsec.out files and writes to toten.dat. After that, the toten.dat will be read.  
`iz` represents the index of layers to be calculated. For the simple scenario of computing 3 z values for the tip, `iz` should be set to 2 as default. Add the `atom` option to display the atom positions. Only the top layer within 1 Å is plotted. Add the `tilt` option to use the tilt correction. Add the `bohr` option to use the Bohr as the length unit.  
