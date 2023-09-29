# UWLCM_plotting
a set of scripts for plotting output of the UWLCM model

CONTENTS

plotting:

- /drawbicyc - C++ program for processing UWLCM output; can calculate and plot series, profiles, variable snapshots and spectra; also has C++ programs for averaging and comparing series or profiles from multiple runs, although this is probably better done in Python

- /cases - Python scripts for comparing series and profiles calculated using drawbicyc; specific to each LES case

- /Energy_spectrum - Python scripts for plotting and comparing Fourier spectra of variables; works directly on UWLCM output

- /histograms - Python scripts for plotting and comparing histograms with spatial and/or temporal distribution of variables; also capable of plotting scatter plots of correlation between two variables; works directly on UWLCM output

- /NC_vs_AF - Python scripts for plotting scatter plots of number concentration of cloud droplets vs adiabatic fraction, or vs liquid water content; works directly on UWLCM output; NOTE: redundant to scatter plots in /histograms?

- /papers - Python scripts for making figures used in our papers

- /Size_spectra - Python scripts for comparing droplet size distributions; works directly on UWLCM output 


helpers:

- /Matplotlib_common - Python scripts for reading drawbicyc output; also plotting comparison of series and profiles from drawbicyc

- /UWLCM_plotters - C++ classes for reading UWLCM output; used in drawbicyc

