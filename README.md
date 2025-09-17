# DFT scripts

This collection of tools is designed to analyze and visualize data from first-principles electronic structure calculations. It includes features such as plotting band structures, converting atomic structures, analyzing the density of states (DOS), and preparing inputs for calculations. Various commonly used calculation packages are supported.

## Installation
Clone the repository and use it as is. The Python scripts are written in Python 3. For convenience, you may add the [src](src) directory to your PATH. To install the required dependencies, use the provided [requirements.txt](requirements.txt) file:
```bash
pip install -r requirements.txt
```

In newer versions of Linux that enforce [PEP 668](https://peps.python.org/pep-0668/), the system's Python environment is treated as externally managed. This prevents using `pip` to install packages directly into the system environment to avoid conflicts with the package manager. To install Python packages independently, create a virtual environment:
```bash
sudo apt install python3.12-venv
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Alternatively, if you prefer to install the required packages directly into the system's Python environment without using `pip`, you can use your package manager, such as `apt`:
```bash
sudo apt install python3-matplotlib python3-numpy python3-scipy python3-yaml
```

## Documentation

📖 Descriptions and usage instructions for functionalities are available on the [Wiki](https://github.com/tangzhao20/dftscr/wiki).  

## Citation

If you use the [iterative tilt correction method](https://github.com/tangzhao20/dftscr/wiki/AFM-simulation-preparation) in the AFM simulations, please cite:  
* Zhao Tang, Dingxin Fan, and James R. Chelikowsky, *Real space simulation for state-resolved high-resolution atomic force microscopy of defects in monolayer h-BN*, [Physical Review Materials **9**, 086201](https://doi.org/10.1103/ncc2-rhmb) (2025).  
