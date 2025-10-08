**Scripts:** [mag_pol.py](https://github.com/tangzhao20/dftscr/blob/main/src/mag_pol.py), [mae.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae.py), [mae_bands.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae_bands.py)

## Saturation magnetic polarization

Saturation magnetic polarization (*J*<sub>s</sub>) is achieved when a material is fully magnetized, forming a uniform single-domain state throughout the structure. Within DFT, this quantity can be determined from a spin-polarized SCF calculation using the total magnetic moment (*m*<sub>tot</sub>) and the unit cell volume:  

$$
J_\mathrm{s} = \mu_0 \frac{m_\mathrm{tot}}{\mathrm{volume}}
$$

The [mag_pol.py](https://github.com/tangzhao20/dftscr/blob/main/src/mag_pol.py) script can be used to calculate *J*<sub>s</sub>:  

```bash
python3 mag_pol.py
```

In this context, the symbol *J* denotes the magnetic polarization and should not be confused with current density. In some references, *I* is used in place of *J*.  

## Magnetic anisotropy energy

Magnetic anisotropy energy (MAE) measures the energy difference when the magnetization of a material is oriented along different crystallographic axes. For simplicity, this work assumes that the easy, medium, and hard axes each align with one of the *x*, *y*, or *z* directions. Conventionally, the MAE can be calculated by performing three separate noncollinear calculations (using VASP as an example). However, the calculation can be accelerated by using second-order perturbation theory (2PT).

In the [mae.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae.py) code, the SOC energy correction using 2PT along direction *i* is:

$$
E_i = \sum_{k} w_k 
\sum_{\sigma_1,\sigma_2} \sigma_1 \sigma_2 
\sum_{b_1,b_2} f_{\sigma_1 b_1 k} (1 - f_{\sigma_2 b_2 k}) 
\frac{E_{\sigma_2 b_2 k} - E_{\sigma_1 b_1 k}}{(E_{\sigma_2 b_2 k} - E_{\sigma_1 b_1 k})^2 + \eta^2} 
\left|\langle \sigma_1 b_1 k | L_i | \sigma_2 b_2 k \rangle\right|^2
$$

Then, the code determines the easy and medium axes and calculates first-order magnetic anisotropy constant *K*<sub>1</sub>:

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

Currently, the code supports the *d* orbitals of Fe and Co, with their spin-orbit coupling constants taken from the following reference:
* M. Blanco-Rey, J. I. Cerda, and A. Arnau, *Validity of perturbative methods to treat the spin–orbit interaction: application to magnetocrystalline anisotropy*, [New Journal of Physics **21**, 73054](https://dx.doi.org/10.1088/1367-2630/ab3060) (2019).  

Since the MAE formula from 2PT is a summation, it can be decomposed to find the contribution of individual bands. [mae_bands.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae_bands.py) calculates the band-resolved MAE contributions. This helps identify the electronic origin of magnetic anisotropy. The code is currently under development.

To run the code, first perform a band structure calculation with `LORBIT = 14`. Then:
```bash
mae_bands.py med easy (E1) (E2)
```

