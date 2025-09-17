**Script:** [dos.py](https://github.com/tangzhao20/dftscr/blob/main/src/dos.py)  

To plot DOS:  
```bash
python3 dos.py (v) package (E1) (E2)
```
* `v` indicates plot vertically   
* `E1` and `E2` define the energy range, as described in the band structure section

**Supported formats:** VASP, QE  

**VASP input:** DOSCAR  
**QE input:** \*.dos

To plot projected DOS:  
```bash
python3 dos.py (v) package atoms orbitals (E1) (E2)
```
**Support format:** QE (qeproj)  

**QE inputs:** \*.dos \*.xml \*.pdos\_atm#\*(\*)\_wfc#\*(\*)

