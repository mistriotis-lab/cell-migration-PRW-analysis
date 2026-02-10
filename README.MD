# Cell Motility Analysis (PRW Model)

This repository contains MATLAB scripts for calculating Mean Squared Displacement (MSD) and fitting a Persistent Random Walk (PRW) model to 2D cell trajectories.
This workflow was developed for and used in **Akinpelu et al. (2026)**.

## Overview

The primary script process tracks cell data (XY coordinates over time) to extract:
- **Speed (S)**: The rate of cell movement
- **Persistence Time (P)**: The average time a cell maintains its direction
- **α (Alpha)**: The anomalous diffusion exponent (Power Law fit)

## Repository Structure

This repository is organized as follows:
- `/src`: contains the core MATLAB functions
- `/example_data`: One sample `.xlsx` file and one sample `.csv` file demonstrating required input format
- `NOTICE`: Detailed attribution of reused PRW fitting logic

## Usage

1. **Dependencies**: Requires MATLAB and the **Curve Fitting Toolbox**

2. **Input Format**: Data should be formatted as `.csv` or `.xlsx` files with the following columns:
    - `TID`: Integer identifiers to determine which cell track data belongs to
    - `x [micron]`: X coordinates for cells at each analyzed frame
    - `y [micron]`: Y coordinates for cells at each analyzed frame
    - `t [min]` **or** `t [mins]`: Time indices for each frame that was sampled

3. **Running the Analysis**: Processing may be run on a batch of files or a single file. 
    - For a single file, you may either run the file or call `single_file_msd_process.m` from the MATLAB command line.
    You will be prompted by a GUI to select the file to be processed.
    - For multiple files within a single directory, run or call `multi_file_msd_process.m` from the MATLAB command line.
    You will be prompted by a GUI to select a folder containing the data files to be processed. 

    In both cases, a subdirectory titled `temp_out` will be generated containing the output files.
    All output files have names derived from the data file that was analyzed, with the string `PROCESSED_` prepended to the file name.

## Mathematical Implementation & Modifications

To ensure statistical goodness in fits, this implementation uses only the first half of MSD data for fitting.
This is done to minimize the influence of noise at high lag times. 

The data is fit to the PRW model:

$$MSD(\tau) = nS^2P[\tau-P(1-exp(-\tau/P))] + 4\sigma ^2$$

## Attributions & Citations

This work builds upon the PRW framework established by the Wirtz lab.
If you use this code, please cite:
1. AKinpelu et al. 2026
2. P.-H. Wu, A. Giri, and D. Wirtz, “Statistical analysis of cell migration in 3D using the anisotropic persistent random walk model,” *Nat Protoc*, vol. 10, no. 3, pp. 517–527, Mar. 2015, doi: [10.1038/nprot.2015.030](10.1038/nprot.2015.030).
3. P.-H. Wu, J. M. Phillip, W. Du, A. Forjaz, P. R. Nair, and D. Wirtz, “Methods to analyze cell migration data: fundamentals and practical guidelines,” *Nat Methods*, vol. 23, no. 1, pp. 43–55, Aug. 2025, doi: [10.1038/s41592-025-02935-5](10.1038/s41592-025-02935-5).