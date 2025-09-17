**Scripts:** [afm.py](https://github.com/tangzhao20/dftscr/blob/main/src/afm.py), [afm.sh](https://github.com/tangzhao20/dftscr/blob/main/src/afm.sh), [afmplot.py](https://github.com/tangzhao20/dftscr/blob/main/src/afmplot.py)

[afm.py](https://github.com/tangzhao20/dftscr/blob/main/src/afm.py) and [afm.sh](https://github.com/tangzhao20/dftscr/blob/main/src/afm.sh) prepare the files and job directories for the AFM simulation. The Python code can be used to prepare the structure files by
```bash
python3 afm.py (vasp)
```

**Inputs:** afm.in, tip.xyz, sample.parsec\_st.dat, (\*\_POTRE.DAT)  

**Outputs:** parsec\_st\_\*\_\*.dat, (parsec\_st\_spot.dat,) manual\_\*\_\*.dat, steps.dat  

An example of [afm.in](https://github.com/tangzhao20/dftscr/blob/main/data/inputs/afm.in) can be found in the [data/inputs](https://github.com/tangzhao20/dftscr/tree/main/data/inputs) directory. If `fdet` is set in afm.in, parsec\_st\_spot.dat will be generated, and the pseudopotential files \*\_POTRE.DAT will be read to calculate Nb. Add `spin` in afm.in to enable spin polarization in the simulations. Add `vasp` at the command line to write an example structure (tip + sample) in VASP format.  

The Bash script can be used to create directories and prepare files by
```bash
afm.sh seq 
```

**Inputs:** afm.in, job.sh, parsec.in.head, parsec\_st\_\*\_\*.dat, (parsec\_st\_spot.dat,) manual\_\*\_\*.dat, steps.dat  

**Outputs:** seq\_\*\_\*, (spot,) \*/parsec.in, \*/job.sh, sbatch.log  

If the `fdet` option is set, enter the spot directory and submit the job manually before iteratively submitting all other jobs. If `fdet` is not set, use the following command to submit the jobs:
```bash
afm.sh sbatch
```

After all the calculations are done, use the Python script [afmplot.py](src/afmplot.py) to make the plots:  
```bash
python3 afmplot.py (iz) (atom) (tilt) (bohr) (verbose) (toten)
```

**Inputs:** afm.in, steps.dat, (toten.dat or seq\_\*\_\*/parsec.out)  

**Outputs:** afm\_\*.png, (toten.dat)  

In the first round, this code reads the total energies from seq\_\*\_\*/parsec.out files and writes to toten.dat. After that, the toten.dat will be read. `iz` represents the index of layers to be calculated. For the simple scenario of computing 3 z values for the tip, `iz` should be set to 2 as default. Add the `atom` option to display the atom positions. Only the top layer within 1 Å is plotted. Add the `tilt` option to use the tilt correction. Add the `bohr` option to use the Bohr as the length unit.  

An iterative tilt correction method is implemented in this code, which solves the lateral shift using $\mathit{\Delta}(x)=F(x+\mathit{\Delta}(x))/k$. The maximum number of iterations, `niter`, and the damping factor, `alpha`, can be set in the input file `afm.in`. This method is described in:  
* Zhao Tang, Dingxin Fan, and James R. Chelikowsky, *Real space simulation for state-resolved high-resolution atomic force microscopy of defects in monolayer h-BN*, [Physical Review Materials **9**, 086201](https://doi.org/10.1103/ncc2-rhmb) (2025).  

The conventional spring model, which assumes $\mathit{\Delta}(x)=F(x)/k$, can be reproduced as a single-shot run of the iterative method by setting `niter 1` and `alpha 1` in the input.

