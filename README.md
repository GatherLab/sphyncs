# sphyncs <a href="https://doi.org/10.5281/zenodo.8121135"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.8121135.svg" alt="DOI"></a>
Software for Processing Hyperspectral Confocal Scans of Laser Particles

This software package contains the Python package for processing hyperspectral confocal datasets, 
as well as the MATLAB Asymptotic Expansion script for mathematical fitting of spherical laser particles.
Custom hardware control code is available upon request due to licensing requirements. 

System requirements: The software packages have been tested on Windows 10.x, using MATLAB R2021a or R2023a, and Python 3.9.13.

MATLAB dependencies: Parallel computing toolbox 7.4, global optimization toolbox 4.5, optimization toolbox 9.1

Python dependencies: pandas 1.4.4, numpy 1.21.5, scipy 1.9.1, matplotlib 3.5.2, imageio 2.19.3, napari 0.4.17, pims 0.6.1, trackpy 0.6.1
No additional non-standard hardware is required to runour processing functions.

Installation guide

MATLAB: Download and install MATLAB and enter the licensing information. Open MATLAB and navigate to 'Environments' to install add-ons (dependencies).
Installing MATLAB and all dependencies will take approximately 30 minutes.

Python: We recommend using the Anaconda distribution, which can be downloaded from the Anaconda website.
Addtional packages (dependencies) can be viewed, updated, or installed under 'Environments' in the Anaconda Navigator, or directly from the anaconda prompt.
Installing Python and all dependencies will take approximately 30 minutes.

Demo: Download all code files from the GitHub/Zenodo repository and ensure all Python (.ipynb/.py) files remain in the same directory. 
Download the Demo source data from the Zenodo repository and extract all files to a directory of your choice.
Open the Demo Jupyter Notebooks and the MATLAB script and update the directory in the code to match the location of the files on your PC.
The code can now be run, reproducing the results published in the protocol 'Hyperspectral Confocal Imaging for High-Throughput Readout and Analysis of Bio-Integrated Laser Particles'.
The expected results are also displayed as graphics in the original Demo notebooks. Note that due to file size limits, only the partial data set corresponding to SI Code 3 was uploaded.
Approximate Demo Runtimes: // SI Code 1: 3 min // SI Code 2: < 1 hour (Python) + < 5 hours (MATLAB) // SI Code 3: < 10 min // SI Code 4: < 1 min

Instructions for Use: Replace the directory, data set name, and log file name by references to your own data.
For optimal results, the keyword arguments of the processing and fitting functions may require updating.
Refer to the Supplementary Information of the protocol 'Hyperspectral Confocal Imaging for High-Throughput Readout and Analysis of Bio-Integrated Laser Particles' for more details.
