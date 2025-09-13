**Script:** [posconvert.py](https://github.com/tangzhao20/dftscr/blob/main/src/posconvert.py)  

To convert the format of the atomic structure:  
```bash
python3 posconvert.py package1 package2 (filename1)
```
**Supported formats:** VASP, QE (.in or .xml), Paratec, PARSEC, .xyz  

**VASP input:** POSCAR  
**QE input:** nscf.in/scf.in/relax.in  
**QExml input:** \*.xml  
**Paratec input:** input  
**PARSEC input:** parsec.in  
**.xyz input:** \*.xyz  

**Optional input:** posconvert.in  

If the optional input file posconvert.in exists, this code modifies the structure before output. The supported operations include move, rotate, flip, vacuum, supercell, add atom, or delete atom. An example of [posconvert.in](https://github.com/tangzhao20/dftscr/blob/main/data/inputs/posconvert.in) can be found in [data/inputs](https://github.com/tangzhao20/dftscr/tree/main/data/inputs) directory.

