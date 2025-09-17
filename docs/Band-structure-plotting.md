**Script:** [bands.py](https://github.com/tangzhao20/dftscr/blob/main/src/bands.py)  

For basic band structure plotting, use:  
```bash
python3 bands.py package (E1) (E2)
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
python3 bands.py package atoms orbitals (E1) (E2)
```
**Supported formats:** VASP (vaspproj), QE (qeproj)  

**VASP inputs:** EIGENVAL, KPOINTS, POSCAR, PROCAR, (DOSCAR)  
**QE inputs:** \*.xml, kpath.in, projwfc.out

To plot band structures from Wannier90:  
```bash
python3 bands.py wan (pad=*) (E1) (E2)
```
**Inputs:** nscf.in, ../bands/\*.xml, \*\_band.kpt, \*\_band.dat, kpath.in

