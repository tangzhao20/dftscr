This directory contains examples of relaxed atomic structures.  

QE is used for the relaxation processes, with the [optB86b-vdW](https://doi.org/10.1103/PhysRevB.83.195131) functional applied in the calculations. The `nc-sr-05_pbe_stringent_upf` version of pseudopotentials from [PseudoDojo](http://www.pseudo-dojo.org/) is utilized, with the normal level of cutoff hints applied. The convergence threshold for self-consistency is set to `conv_thr = 1.d-9`. The k-points are selected to ensure that `Nk_i * a_i > 30 Å`.  

The structures here include [monolayer h-BN](hbn.vasp) and [bulk silicon](silicon.vasp).  

