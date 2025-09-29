**Scripts:** [mae.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae.py), [mae_bands.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae_bands.py)

Magnetic anisotropy energy (MAE) measures the energy difference when the magnetization of a material is oriented along different crystallographic axes. For simplicity, this work assumes that the easy, medium, and hard axes each align with one of the *x*, *y*, or *z* directions. Conventionally, the MAE can be calculated by performing three separate noncollinear calculations (using VASP as an example). However, the calculation can be accelerated by using second-order perturbation theory (2PT).

In the [mae.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae.py) code, the SOC energy correction using 2PT along direction *i* is:

$$
E_i = \sum_{k} w_k 
\sum_{\sigma_1,\sigma_2} \sigma_1 \sigma_2 
\sum_{b_1,b_2} f_{\sigma_1 b_1 k} (1 - f_{\sigma_2 b_2 k}) 
\frac{E_{\sigma_2 b_2 k} - E_{\sigma_1 b_1 k}}{(E_{\sigma_2 b_2 k} - E_{\sigma_1 b_1 k})^2 + \eta^2} 
\left|\langle \sigma_1 b_1 k | L_i | \sigma_2 b_2 k \rangle\right|^2
$$

Then the code determines the easy and medium axes, and calculates $K_1$:

$$
K_1 = \frac{E_{\mathrm{med}}-E_{\mathrm{easy}}}{\mathrm{volume}}
$$

To calculate MAE, first run a VASP collinear SCF calculation with the following settings:
- `LORBIT = 14` to output orbital projections with phase.  
- Initialize `MAGMOM` properly (recommended: 4 for Fe and Co).  
- The current version does not support symmetry operations with rotations higher than threefold. If the atomic structure contains such symmetries, set `ISYM = 0` to avoid using them.  

Then run:
```bash
python3 mae.py
```

Currently, the code support the *d* orbitals of Fe and Co.

---

Since the MAE formula from 2PT is a summation, it can be decomposed to find the contribution of individual bands. [mae_bands.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae_bands.py) calculates the band-resolved MAE contributions. This helps identify the electronic origin of magnetic anisotropy. The code is currently under development.

To run the code, first perform a band structure calculation with `LORBIT = 14`. Then:
```bash
mae_bands.py med easy (E1) (E2)
```

