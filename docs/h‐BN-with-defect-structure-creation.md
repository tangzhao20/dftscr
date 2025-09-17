**Script:** [hbndef.py](https://github.com/tangzhao20/dftscr/blob/main/src/hbndef.py)

This script generates h-BN monolayer or flake structures with or without a defect.

Generate 2D h-BN structures of a hexagonal cell:
```bash
python3 hbndef.py N (defect)
```

Generate 2D h-BN structures of a rectangular cell:
```bash
python3 hbndef.py rec r (defect)
```

Generate h-BN flake structures:
```bash
python3 hbndef.py flake r (defect)
```

* `N` defines the size of the supercell.  
* `r` denotes either the radius of the h-BN flake or half the minimum defect–defect distance.  
* `(defect)` is optional. If provided, it indicates the defect name, such as CBVN. Make sure to use the correct capitalization for atomic symbols.  

