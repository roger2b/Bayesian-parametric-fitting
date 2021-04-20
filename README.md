# Bayesian-parametric-fitting
Fast and accurate extraction of cosmic information signal (eg. from the CMB)

## codes
foreground_cleaning_pixel.py - main code to read in uncleaned CMB data and noise covariance matrices and perform a Bayesian parametric fitting. 
plot_cleaned_map.py          - code to read outpouts from foreground_cleaning_pixel.py to plot results and output maps (all as .txt files) 

## nest steps
This repo is under development and paper is work in progress. 
- Read flexible number of freq channels
- Read constant or full noise covariance matrices
- choose to include pixel-pixel correlations (for now noise is assumed to be diagonal) 
- parallelisation of computations

## citation
Please cite de Belsunce et al. (2021) in prep. if you use this or parts of this code
