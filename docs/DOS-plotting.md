**Script:** [dos.py](https://github.com/tangzhao20/dftscr/blob/main/src/dos.py)  

To plot DOS:  
```bash
python3 dos.py (v) (o) package (E1) (E2)
```

To plot projected DOS (PDOS):  
```bash
python3 dos.py (v) (o) package atom_list orb_list (E1) (E2)
```

**Options**  
* `v`: plot vertically (e.g., for side-by-side band structure comparison).  
* `o`: output data to .dat files.  
* `E1` `E2`: energy range. Default => -5~5; `E1` => -E1~E1; `E1` `E2` => E1~E2.  

**Supported formats:** VASP, QE  

| Code | Mode | `package` arg | Required files |
| :--- | :--- | :--- | :--- |
| **VASP** | DOS | `vasp` | `DOSCAR` |
| | PDOS | `vaspproj` | `DOSCAR`, `POSCAR`, `PROCAR` |
| **QE** | DOS | `qe` | `*.dos`, `*.xml` |
| | PDOS | `qeproj` | `*.dos`, `*.xml`, `*.pdos_atm#*(*)_wfc#*(*)` |
