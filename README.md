# Bayesian-parametric-fitting
Fast and accurate component seperation of primary CMB and Galactic polarised foregrounds

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
Please cite de Belsunce et al. (2021) if you use this or parts of this code

@ARTICLE{2022arXiv220513968D,
       author = {{de Belsunce}, Roger and {Gratton}, Steven and {Efstathiou}, George},
        title = "{Testing for spectral index variations in polarised CMB foregrounds}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Cosmology and Nongalactic Astrophysics},
         year = 2022,
        month = may,
          eid = {arXiv:2205.13968},
        pages = {arXiv:2205.13968},
archivePrefix = {arXiv},
       eprint = {2205.13968},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022arXiv220513968D},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

