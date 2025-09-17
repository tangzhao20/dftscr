**Scripts:** [kconvert.py](https://github.com/tangzhao20/dftscr/blob/main/src/kconvert.py)

To convert the format of the *k*-point path:  
```bash
python3 kconvert.py package1 package2 (N)
```
**Supported formats:** VASP, QE, kpath.in

**VASP input:** KPOINTS  

Examples of the kpath.in file can be found in the [data/kpaths](https://github.com/tangzhao20/dftscr/tree/main/data/kpaths) directory.

A common use of this code is to create the *k*-point path for QE band structure calculations:  
```bash
python3 kconvert.py kpathin qe (N)
```

