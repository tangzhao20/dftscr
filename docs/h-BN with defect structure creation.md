**Script:** [hbndef.py](https://github.com/tangzhao20/dftscr/blob/main/src/hbndef.py)

This script generates h-BN monolayer or flake structures with or without a defect. To create a 2D h-BN structure, use the following command:
```bash
python3 hbndef.py N (defect)
```
For an h-BN flake, use:
```bash
python3 hbndef.py flake r (defect)
```
* `N` is the supercell size.  
* `r` is the radius of the h-BN flake.  
* `defect` specifies the defect, such as CBVN. Ensure atomic symbols are correctly capitalized.  

