# GaMDWE - Gaussian Accelerated Molecular Dynamics - Weighted Ensemble (GaMDWE) 

## Installation and Setup Instructions :
* Make sure [anaconda3](https://www.anaconda.com/) is installed on the local machine. Go to the  [download](https://www.anaconda.com/products/individual) page of anaconda3 and install the latest version of anaconda3.
* Create a new conda environment with python = 3.8 and install the package with the following commands in the terminal: 
```bash
conda create -n gamdwe python=3.6
conda activate gamdwe
conda install -c conda-forge curl
conda install -c conda-forge matplotlib
conda install -c conda-forge openmm
conda install -c conda-forge seaborn
conda install -c conda-forge pandas 
conda install -c ambermd pytraj
conda install -c omnia ambertools
```
