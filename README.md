# PhageCounting
This repository contains the code and data required to reproduce the results in the manuscript “Using bacterial population dynamics to count phages and their lysogens”.

# 1. Folder Structure
Each folder is named according to the figure(s) it contains. Inside each folder, you will find the following subfolders:

* `data` (optional): Contains the raw data used to generate the plots.
* `script`: Contains the Python/MATLAB scripts used to generate the plots.
* `tmp` (optional): Stores intermediate files generated during code execution, such as ensembles of model parameters.

# 2. Running the Codes
Each script is named according to the figure panel(s) it corresponds to. The scripts are either MATLAB files (for Fig. 4, Figs. 5b-e, and Sfig. 17) or Jupyter notebooks (for other figures).

## 2.1. Jupyter Notebooks
The scripts were developed using Python 3.8.5. To run the Jupyter notebooks, you need the following packages:

* numpy
* scipy
* matplotlib
* pandas
* emcee
* lmfit

To store the output figures, you need to create a `output` subfolder, inside the main folder of the figure. 

For scripts involving MCMC (Fig. 3c, Fig. 5f, Sfig. 7, Sfig. 12, Sfig. 17, Sfig. 19, and Sfig. 20), we have provided the sampled parameters from previous runs in `.csv` in the `tmp` folder, which can be directly loaded into the notebooks using `pd.read_csv`.

## 2.2. MATLAB
The scripts were developed using MATLAB versions 2020a–2023b (64-bit). To run the MATLAB scripts:

* Set the folder containing the `.m` files as the working directory for the current MATLAB session.
* Add the relevant folders and subfolders to the search path of the current MATLAB session.
* Modify the paths in the MATLAB scripts to correctly load the data files.
* The sections of the code must be run from top to bottom (in the order in the script file). 

# 3. Reference
Please cite the following paper:
Geng, Yuncong, et al. "Using population dynamics to count bacteriophages and their lysogens." *bioRxiv* (2023): 2023-10.

# 4. Contact Information
- Yuncong Geng (yuncong5@illinois.edu)
- Thu Vu Phuc Nguyen (thuvpnguyen@princeton.edu)
- Ido Golding (igolding@illinois.edu)
