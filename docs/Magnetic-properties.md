**Scripts:** [mag_pol.py](https://github.com/tangzhao20/dftscr/blob/main/src/mag_pol.py), [mae.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae.py), [mae_bands.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae_bands.py), [mae_density.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae_density.py)

## Saturation magnetic polarization

Saturation magnetic polarization (*J*<sub>s</sub>) is achieved when a material is fully magnetized, forming a uniform single-domain state throughout the structure. Within DFT, this quantity can be determined from a spin-polarized SCF calculation using the total magnetic moment (*m*<sub>tot</sub>) and the unit cell volume:  

$$
J_\mathrm{s} = \mu_0 \frac{m_\mathrm{tot}}{\mathrm{volume}}
$$

The [mag_pol.py](https://github.com/tangzhao20/dftscr/blob/main/src/mag_pol.py) script can be used to calculate *J*<sub>s</sub>:  

```bash
mag_pol.py
```

In this context, the symbol *J* denotes the magnetic polarization and should not be confused with current density. In some references, *I* is used in place of *J*.  

## Magnetic anisotropy energy

Magnetic anisotropy energy (MAE) measures the energy difference when the magnetization of a material is oriented along different crystallographic axes. For simplicity, this work assumes that the easy, medium, and hard axes each align with one of the *x*, *y*, or *z* directions. Conventionally, the MAE can be calculated by performing three separate noncollinear DFT calculations (using VASP as an example). These calculations can be sped up by using second-order perturbation theory (2PT) as implemented in [mae.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae.py).  

The 2PT expression contains a term involving the inverse of the energy differences between band pairs:

$$
\mathrm{S}(\varepsilon_1,\varepsilon_2) = \frac{1}{\varepsilon_1-\varepsilon_2}
$$

This term diverges at $\varepsilon_1=\varepsilon_2$, in which case the system should be treated with degenerate perturbation theory. The following discussion assumes the degenerate effect to be negligible upon integration over the Brillouin zone. A Lorentzian-type smearing is applied as a replacement to avoid the divergence:  

$$
\mathrm{S}_\eta(\varepsilon_1,\varepsilon_2) = \frac{\varepsilon_1-\varepsilon_2}{(\varepsilon_1-\varepsilon_2)^2+\eta^2}
$$


An angular momentum term is defined as:  

$$
L_{i \sigma_1 b_1 \sigma_2 b_2 k} = \sum_a \xi_a \langle \sigma_1 b_1 a k | L_i | \sigma_2 b_2 a k \rangle
$$
 
From this, the SOC energy correction along direction *i* is given by:

$$
E_i = \sum_{k} w_k 
\sum_{\sigma_1,\sigma_2} \sigma_1 \sigma_2 
\sum_{b_1,b_2} f_{\sigma_1 b_1 k} (1 - f_{\sigma_2 b_2 k}) 
\mathrm{S}_\eta(\varepsilon_{\sigma_2 b_2 k}, \varepsilon_{\sigma_1 b_1 k})
\left|L_{i \sigma_1 b_1 \sigma_2 b_2 k}\right|^2
$$

Based on these directional energy calculations, the code identifies the easy and medium axes to obtain the MAE:

$$
\mathrm{MAE} = E_{\mathrm{med}}-E_{\mathrm{easy}}
$$

and the first-order magnetic anisotropy constant *K*<sub>1</sub>:

$$
K_1 = \frac{\mathrm{MAE}}{\mathrm{volume}}
$$

To calculate MAE, first perform a collinear SCF calculation using VASP with the following settings:
- `LORBIT = 14` to output orbital projections with phase.  
- Initialize `MAGMOM` properly (recommended: 4 for Fe and Co).  
- The current version does not support symmetry operations with rotations higher than threefold. If the atomic structure contains such symmetries, set `ISYM = 0` to avoid using them.  
- High `NBANDS` may cause incorrect results. The reason is under investigation.  

Execute the following command to perform the calculations for MAE:
```bash
mae.py [Nb=<int>]
```
**Nb**: Number of bands to include. Default: all available bands.  

Currently, the code supports the *d* orbitals of Fe and Co, and their spin-orbit coupling constants (*ξ<sub>a</sub>*) are taken from the following reference:
* M. Blanco-Rey, J. I. Cerda, and A. Arnau, *Validity of perturbative methods to treat the spin–orbit interaction: application to magnetocrystalline anisotropy*, [New Journal of Physics **21**, 73054](https://dx.doi.org/10.1088/1367-2630/ab3060) (2019).  

## Decomposition of MAE

Since the MAE from 2PT takes the form of a summation, it can be analyzed by decomposing its contributions. This allows one to identify the features on the band structure that dominate the anisotropy. The following decomposition approaches are implemented in the scripts.  

### *k*- and band-resolved decomposition
The [mae_bands.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae_bands.py) script allows *k*- and band-resolved decomposition of MAE and produces a plot in the style of a projected band structure.  Essentially, this quantity follows the MAE expression but includes only a single band summation at a fixed *k*-point.  

For spin *σ*<sub>1</sub> and band *b*<sub>1</sub> at *k*:  

$$
\mathrm{MAE}_{\sigma_1 b_1 k} = \frac{1}{2} \sum_{\sigma_2} \sigma_1 \sigma_2 \sum_{b_2}
   \left( f_{b_1 \sigma_1 k} - f_{b_2 \sigma_2 k} \right) 
   \mathrm{S}_\eta(\varepsilon_{\sigma_2 b_2 k}, \varepsilon_{\sigma_1 b_1 k})
   \left( \left| L_{i_\mathrm{med} \sigma_1 b_1 \sigma_2 b_2 k} \right|^2 -
   \left| L_{i_\mathrm{easy} \sigma_1 b_1 \sigma_2 b_2 k} \right|^2 \right)
$$

where *i*<sub>med</sub> and *i*<sub>easy</sub> refer to the medium and easy axis directions.

To run the code, first perform a band structure calculation with `LORBIT = 14`. Then execute:
```bash
mae_bands.py med easy (E1) (E2)
```

`med` and `easy` should be chosen from `x`, `y`, or `z`.  

This decomposition method provides the most detailed information. The script presents positive and negative MAE contributions as red and green scatter points, plotted on top of the band structure, where the two spin channels are shown as dark blue and orange lines. The output is saved as `mae_bs.png`. However, the combination of four colors and overlapping features makes the plot overly complex, so further summarization may be needed for clarity.  

### *k*-resolved decomposition

The *k*-resolved decomposition sums contributions over bands to show the MAE from each *k*-point along the band-structure path. This highlights which regions of the Brillouin zone dominate the anisotropy.  

$$
\mathrm{MAE}_{\sigma_1 k} = \sum_{b_1} \mathrm{MAE}_{\sigma_1 b_1 k}
$$

A plot of this decomposition is generated by the [mae_bands.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae_bands.py) script, alongside the previously described [*k*- and band-resolved decomposition plot](#k--and-band-resolved-decomposition). No additional step is needed. The output is saved as `mae_bs_sum.png`, with MAE contributions shown in arbitrary units.  

### Energy-resolved MAE density

To provide a complementary perspective, [mae_density.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae_density.py) generates a DOS-style distribution of MAE as a function of energy:

$$
\mathrm{MAE}_{\sigma_1}(E) = \sum_{b_1 k} w_k \mathrm{MAE}_{\sigma_1 b_1 k} 
\mathrm{G}_\mathit{\Delta}(E-\varepsilon_{\sigma_1 b_1 k})
$$

where $\mathrm{G}_\mathit{\Delta}$ is the Gaussian smearing function:

$$
\mathrm{G}_\mathit{\Delta}(E-\varepsilon_1) = \frac{1}{\mathit{\Delta}\sqrt{2\pi}} 
\mathrm{exp} \left( -\frac{(E-\varepsilon_1)^2}{2\mathit{\Delta}^2} \right)
$$

The total MAE can be recovered by integrating $\mathrm{MAE}_{\sigma_1}(E)$ over energy within a sufficiently large window and summing over spins.

This script should be run after the same DFT calculations for [mae.py](https://github.com/tangzhao20/dftscr/blob/main/src/mae.py), which are performed on *k*-grids with `LORBIT = 14`:  
```bash
mae_density.py med easy (E1) (E2)
```

See the [MAE section](#magnetic-anisotropy-energy) for details on the requirements. The output is saved as `mae_density.png`, with MAE density shown in arbitrary units.  
